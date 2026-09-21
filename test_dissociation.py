"""
test_dissociation.py — the outside of dissociation.py.

unittest, because the module Depends on nothing but the standard library
and a test suite is no reason to invite a guest.

Every expectation is Spelled: a textbook number, a chemical invariant
written out here, or a limit the chemistry already knows.
"""

import dataclasses
import math
import unittest

import dissociation as d


# --------------------------------------------------------------------------
# the closed form
# --------------------------------------------------------------------------

class TheMonoproticFormulaSolvesItsOwnQuadratic(unittest.TestCase):

    def test_the_answer_satisfies_the_definition_of_ka(self):
        for ka, c in ((1.8e-5, 0.1), (6.2e-10, 0.05), (1.5e-2, 1.0), (1.0e-19, 0.1)):
            with self.subTest(ka=ka, c=c):
                x = d.monoprotic_h(ka, c)
                self.assertAlmostEqual(x * x / (c - x), ka, delta=ka * 1e-9)

    def test_the_rationalized_form_agrees_with_the_textbook_one(self):
        # (Ka/2)*(sqrt(1 + 4C/Ka) - 1), written here and not read from there.
        for ka, c in ((1.8e-5, 0.1), (6.8e-4, 0.25), (3.0e-8, 0.01)):
            with self.subTest(ka=ka, c=c):
                textbook = (ka / 2.0) * (math.sqrt(1.0 + 4.0 * c / ka) - 1.0)
                self.assertAlmostEqual(d.monoprotic_h(ka, c), textbook, delta=textbook * 1e-9)

    def test_a_strong_acid_dissociates_almost_completely(self):
        # Ka = 1e6 >> C: every proton is free, so [H+] -> C.
        self.assertAlmostEqual(d.monoprotic_h(1e6, 0.01), 0.01, delta=1e-7)

    def test_nothing_dissolved_releases_nothing(self):
        self.assertEqual(d.monoprotic_h(1.8e-5, 0.0), 0.0)

    def test_a_constant_must_be_positive_and_a_concentration_must_not_be_negative(self):
        with self.assertRaises(ValueError):
            d.monoprotic_h(0.0, 0.1)
        with self.assertRaises(ValueError):
            d.monoprotic_h(-1.8e-5, 0.1)
        with self.assertRaises(ValueError):
            d.monoprotic_h(1.8e-5, -0.1)


# --------------------------------------------------------------------------
# water, where it matters
# --------------------------------------------------------------------------

class WaterCompetesWhenTheAcidIsWeakOrDilute(unittest.TestCase):

    def test_the_balance_the_solver_closes_is_the_one_it_promises(self):
        # h = sum_j j*alpha_j*C + Kw/h, checked from outside.
        kas = [7.5e-3, 6.2e-8, 4.2e-13]
        c = 0.1
        h = d.polyprotic_h(kas, c)
        alphas = d.fractions(kas, h)
        released = c * sum(j * a for j, a in enumerate(alphas))
        self.assertAlmostEqual(h, released + d.KW / h, delta=h * 1e-9)

    def test_pure_water_is_neutral(self):
        self.assertAlmostEqual(d.acid_ph(1.8e-5, 0.0, exact=True), 7.0, places=9)

    def test_an_acid_can_never_be_basic_however_dilute(self):
        # 1e-8 M of a strong acid: the naive formula says pH 8, which is absurd.
        naive = d.acid_ph(1e6, 1e-8)
        exact = d.acid_ph(1e6, 1e-8, exact=True)
        self.assertGreater(naive, 7.0)          # the known failure of ignoring water
        self.assertLess(exact, 7.0)
        self.assertGreater(exact, 6.9)

    def test_a_very_weak_acid_at_low_concentration_lands_just_under_neutral(self):
        # HCN 1e-6 M: water carries most of the protons.
        self.assertAlmostEqual(d.acid_ph(6.2e-10, 1e-6, exact=True), 6.987, places=3)
        self.assertGreater(d.acid_ph(6.2e-10, 1e-6), 7.0)   # ignoring water overshoots

    def test_an_absurdly_small_ka_leaves_the_solution_neutral(self):
        self.assertAlmostEqual(d.acid_ph(1.0e-30, 0.1, exact=True), 7.0, places=6)

    def test_water_is_irrelevant_for_a_concentrated_weak_acid(self):
        loose = d.acid_ph(1.8e-5, 0.1)
        tight = d.acid_ph(1.8e-5, 0.1, exact=True)
        self.assertAlmostEqual(loose, tight, places=6)

    def test_acetic_acid_at_a_tenth_molar_is_the_number_in_the_book(self):
        self.assertAlmostEqual(d.acid_ph(1.8e-5, 0.1), 2.875, places=3)

    def test_ammonia_at_a_tenth_molar_is_the_number_in_the_book(self):
        self.assertAlmostEqual(d.base_ph(1.8e-5, 0.1), 11.125, places=3)

    def test_an_acid_and_its_mirror_base_add_up_to_pkw(self):
        self.assertAlmostEqual(
            d.acid_ph(1.8e-5, 0.1) + d.base_ph(1.8e-5, 0.1), 14.0, places=9)

    def test_a_different_pkw_moves_the_neutral_point(self):
        # Kw = 1e-13 (about 60 C): neutral sits at 6.5.
        self.assertAlmostEqual(d.acid_ph(1.8e-5, 0.0, kw=1e-13, exact=True), 6.5, places=9)

    def test_a_different_pkw_moves_the_base_side_too(self):
        # At Kw = 1e-13 a vanishing base leaves neutral water at 6.5, not 7.
        self.assertAlmostEqual(d.base_ph(1.8e-5, 0.0, kw=1e-13, exact=True), 6.5, places=9)
        # and ammonia measured against the same pKw sits one unit lower.
        self.assertAlmostEqual(
            d.base_ph(1.8e-5, 0.1, kw=1e-13), d.base_ph(1.8e-5, 0.1) - 1.0, places=9)

    def test_the_balance_needs_at_least_one_constant_and_a_real_concentration(self):
        with self.assertRaises(ValueError):
            d.polyprotic_h([], 0.1)
        with self.assertRaises(ValueError):
            d.polyprotic_h([1.8e-5], -0.1)
        with self.assertRaises(ValueError):
            d.polyprotic_h([1.8e-5, 0.0], 0.1)


# --------------------------------------------------------------------------
# speciation
# --------------------------------------------------------------------------

class TheFractionsDivideTheWholeSolute(unittest.TestCase):

    def test_they_always_sum_to_one(self):
        kas = [7.4e-4, 1.7e-5, 4.0e-7]
        for exponent in range(0, 15):
            h = 10.0 ** -exponent
            with self.subTest(ph=exponent):
                self.assertAlmostEqual(sum(d.fractions(kas, h)), 1.0, places=12)
                self.assertTrue(all(0.0 <= a <= 1.0 for a in d.fractions(kas, h)))

    def test_one_constant_gives_two_species(self):
        self.assertEqual(len(d.fractions([1.8e-5], 1e-4)), 2)

    def test_at_ph_equal_to_pka_the_pair_is_split_in_half(self):
        acid, base = d.fractions([1.8e-5], 1.8e-5)
        self.assertAlmostEqual(acid, 0.5, places=12)
        self.assertAlmostEqual(base, 0.5, places=12)

    def test_strong_acidity_leaves_only_the_protonated_form(self):
        # h well above Ka1: H3PO4 keeps all three protons.
        alphas = d.fractions([7.5e-3, 6.2e-8, 4.2e-13], 1e4)
        self.assertAlmostEqual(alphas[0], 1.0, places=5)

    def test_strong_alkalinity_leaves_only_the_bare_anion(self):
        # h well below Ka3: only PO4^3- survives.
        alphas = d.fractions([7.5e-3, 6.2e-8, 4.2e-13], 1e-18)
        self.assertAlmostEqual(alphas[-1], 1.0, places=5)

    def test_phosphate_at_neutral_ph_is_mostly_dihydrogen_phosphate(self):
        alphas = d.fractions([7.5e-3, 6.2e-8, 4.2e-13], 1e-7)
        self.assertEqual(alphas.index(max(alphas)), 1)
        self.assertAlmostEqual(alphas[1], 0.6173, places=4)
        self.assertAlmostEqual(alphas[2], 0.3827, places=4)


# --------------------------------------------------------------------------
# buffers
# --------------------------------------------------------------------------

class TheBufferAgreesWithHendersonHasselbalchWhereItHolds(unittest.TestCase):

    def test_an_equimolar_pair_sits_at_the_pka(self):
        self.assertAlmostEqual(d.henderson_hasselbalch(1.8e-5, 0.1, 0.1), d.pk(1.8e-5), places=12)
        self.assertAlmostEqual(d.buffer_ph(1.8e-5, 0.1, 0.1), d.pk(1.8e-5), places=3)

    def test_ten_to_one_moves_the_ph_one_unit(self):
        self.assertAlmostEqual(
            d.henderson_hasselbalch(1.8e-5, 0.01, 0.1) - d.pk(1.8e-5), 1.0, places=12)
        self.assertAlmostEqual(
            d.henderson_hasselbalch(1.8e-5, 0.1, 0.01) - d.pk(1.8e-5), -1.0, places=12)

    def test_the_two_answers_meet_in_the_concentrated_regime(self):
        for acid_c, base_c in ((0.1, 0.1), (0.2, 0.05), (0.5, 0.5), (0.05, 0.2)):
            with self.subTest(acid=acid_c, base=base_c):
                self.assertAlmostEqual(
                    d.buffer_ph(1.8e-5, acid_c, base_c),
                    d.henderson_hasselbalch(1.8e-5, acid_c, base_c),
                    delta=0.01)

    def test_the_two_answers_part_ways_once_the_pair_is_dilute(self):
        # 1e-5 M of each: both members are consumed, and the approximation lies.
        exact = d.buffer_ph(1.8e-5, 1e-5, 1e-5)
        approximate = d.henderson_hasselbalch(1.8e-5, 1e-5, 1e-5)
        self.assertGreater(exact - approximate, 0.4)

    def test_a_buffer_needs_both_members_of_the_pair(self):
        with self.assertRaises(ValueError):
            d.henderson_hasselbalch(1.8e-5, 0.0, 0.1)
        with self.assertRaises(ValueError):
            d.henderson_hasselbalch(1.8e-5, 0.1, 0.0)

    def test_the_exact_buffer_refuses_impossible_numbers(self):
        with self.assertRaises(ValueError):
            d.buffer_h(0.0, 0.1, 0.1)
        with self.assertRaises(ValueError):
            d.buffer_h(1.8e-5, -0.1, 0.1)
        with self.assertRaises(ValueError):
            d.buffer_h(1.8e-5, 0.1, -0.1)

    def test_a_buffer_with_no_base_is_just_the_free_acid(self):
        self.assertAlmostEqual(
            d.buffer_ph(1.8e-5, 0.1, 0.0), d.acid_ph(1.8e-5, 0.1, exact=True), places=6)


# --------------------------------------------------------------------------
# titration
# --------------------------------------------------------------------------

class TheTitrationPassesThroughItsLandmarks(unittest.TestCase):

    ACETIC = 1.8e-5

    def test_the_equivalence_volume_is_a_ratio_of_equivalents(self):
        self.assertAlmostEqual(d.equivalence_volume(0.1, 50.0, 0.1), 50.0, places=12)
        self.assertAlmostEqual(d.equivalence_volume(0.1, 50.0, 0.2), 25.0, places=12)

    def test_a_diprotic_analyte_needs_twice_the_titrant(self):
        self.assertAlmostEqual(d.equivalence_volume(0.1, 50.0, 0.1, protons=2), 100.0, places=12)

    def test_a_titrant_without_concentration_is_refused(self):
        with self.assertRaises(ValueError):
            d.equivalence_volume(0.1, 50.0, 0.0)

    def test_before_any_titrant_the_ph_is_the_free_acid(self):
        self.assertAlmostEqual(
            d.titration_ph([self.ACETIC], 0.1, 50.0, 0.1, 0.0),
            d.acid_ph(self.ACETIC, 0.1, exact=True), places=9)

    def test_at_half_equivalence_the_ph_is_the_pka(self):
        self.assertAlmostEqual(
            d.titration_ph([self.ACETIC], 0.1, 50.0, 0.1, 25.0), d.pk(self.ACETIC), places=2)

    def test_the_equivalence_point_of_a_weak_acid_is_basic(self):
        at_equivalence = d.titration_ph([self.ACETIC], 0.1, 50.0, 0.1, 50.0)
        self.assertGreater(at_equivalence, 8.0)
        self.assertLess(at_equivalence, 9.0)

    def test_excess_strong_base_governs_the_tail(self):
        # 10 mL of 0.1 M NaOH beyond equivalence in 110 mL: [OH-] = 1/110 M.
        expected = 14.0 + math.log10(0.1 * 10.0 / 110.0)
        self.assertAlmostEqual(
            d.titration_ph([self.ACETIC], 0.1, 50.0, 0.1, 60.0), expected, places=3)

    def test_a_total_volume_of_zero_has_no_ph(self):
        with self.assertRaises(ValueError):
            d.titration_ph([self.ACETIC], 0.1, 0.0, 0.1, 0.0)

    def test_the_curve_starts_at_zero_and_ends_where_it_was_told(self):
        curve = d.titration_curve([self.ACETIC], 0.1, 50.0, 0.1, up_to=100.0, points=101)
        self.assertEqual(len(curve), 101)
        self.assertAlmostEqual(curve[0][0], 0.0, places=12)
        self.assertAlmostEqual(curve[-1][0], 100.0, places=12)

    def test_the_curve_walks_through_the_equivalence_point(self):
        curve = d.titration_curve([self.ACETIC], 0.1, 50.0, 0.1, points=101)
        volumes = [v for v, _ in curve]
        self.assertAlmostEqual(volumes[50], 50.0, places=9)
        self.assertAlmostEqual(
            curve[50][1], d.titration_ph([self.ACETIC], 0.1, 50.0, 0.1, 50.0), places=9)

    def test_the_curve_only_rises(self):
        curve = d.titration_curve([self.ACETIC], 0.1, 50.0, 0.1, points=51)
        for (_, before), (_, after) in zip(curve, curve[1:]):
            self.assertGreater(after, before)

    def test_the_jump_is_steepest_at_the_equivalence_point(self):
        curve = d.titration_curve([self.ACETIC], 0.1, 50.0, 0.1, points=101)
        jumps = [curve[i + 1][1] - curve[i][1] for i in range(len(curve) - 1)]
        self.assertEqual(jumps.index(max(jumps)), 49)

    def test_a_diprotic_curve_reaches_past_both_protons(self):
        # The default range is twice the volume both protons need.
        curve = d.titration_curve([5.9e-2, 6.4e-5], 0.1, 50.0, 0.1, points=201)
        self.assertAlmostEqual(curve[-1][0], 200.0, places=9)
        jumps = [(curve[i + 1][1] - curve[i][1], curve[i][0]) for i in range(len(curve) - 1)]
        _, steepest_volume = max(jumps)
        self.assertAlmostEqual(steepest_volume, 99.0, delta=1.0)  # second equivalence

    def test_a_curve_needs_at_least_two_points(self):
        with self.assertRaises(ValueError):
            d.titration_curve([self.ACETIC], 0.1, 50.0, 0.1, points=1)


# --------------------------------------------------------------------------
# constants and the table
# --------------------------------------------------------------------------

class ThePkIsOnlyALogarithm(unittest.TestCase):

    def test_acetic_acid_has_a_pka_near_four_and_three_quarters(self):
        self.assertAlmostEqual(d.pk(1.8e-5), 4.745, places=3)

    def test_the_pair_round_trips(self):
        for value in (1.0e-14, 1.8e-5, 0.5, 1.0, 72.0):
            with self.subTest(value=value):
                self.assertAlmostEqual(d.k_from_pk(d.pk(value)), value, delta=value * 1e-12)

    def test_a_constant_of_zero_or_less_has_no_logarithm(self):
        with self.assertRaises(ValueError):
            d.pk(0.0)
        with self.assertRaises(ValueError):
            d.pk(-1.0)


class TheTableAnswersToEveryNameItWasGiven(unittest.TestCase):

    def test_a_formula_finds_the_species(self):
        self.assertEqual(d.find_species("CH3COOH").name, "acetic acid")

    def test_an_english_name_finds_it_too(self):
        self.assertEqual(d.find_species("acetic acid").formula, "CH3COOH")

    def test_an_accented_spanish_alias_finds_it(self):
        self.assertEqual(d.find_species("ácido acético").formula, "CH3COOH")

    def test_the_accent_is_optional_and_so_is_the_case(self):
        self.assertEqual(d.find_species("  ACIDO ACETICO  ").formula, "CH3COOH")

    def test_a_nickname_finds_it(self):
        self.assertEqual(d.find_species("vinagre").formula, "CH3COOH")
        self.assertEqual(d.find_species("vitamin c").formula, "C6H8O6")

    def test_bases_are_searched_as_well(self):
        self.assertTrue(d.find_species("amoníaco").is_base)

    def test_an_unknown_species_is_refused_by_name(self):
        with self.assertRaises(KeyError):
            d.find_species("unobtainium")

    def test_every_tabulated_constant_is_positive_and_ordered(self):
        for table in (d.ACIDS, d.BASES):
            for formula, s in table.items():
                with self.subTest(formula=formula):
                    self.assertEqual(formula, s.formula)
                    self.assertTrue(all(k > 0 for k in s.k))
                    self.assertEqual(list(s.k), sorted(s.k, reverse=True))

    def test_the_species_reports_its_own_pk(self):
        self.assertEqual(len(d.find_species("H3PO4").pk), 3)
        self.assertAlmostEqual(d.find_species("CH3COOH").pk[0], 4.745, places=3)

    def test_the_repr_names_the_constant_it_carries(self):
        self.assertIn("Ka=", repr(d.find_species("CH3COOH")))
        self.assertIn("Kb=", repr(d.find_species("NH3")))

    def test_a_species_is_frozen(self):
        with self.assertRaises(dataclasses.FrozenInstanceError):
            d.find_species("CH3COOH").name = "something else"


class ThePhEntryPointDispatchesOnTheSpecies(unittest.TestCase):

    def test_an_acid_by_name_is_acidic(self):
        self.assertAlmostEqual(d.ph("CH3COOH", 0.1), 2.875, places=3)

    def test_a_base_by_name_is_basic(self):
        self.assertAlmostEqual(d.ph("NH3", 0.1), 11.125, places=3)

    def test_a_species_object_works_as_well_as_its_name(self):
        acetic = d.find_species("CH3COOH")
        self.assertEqual(d.ph(acetic, 0.1), d.ph("CH3COOH", 0.1))

    def test_a_polyprotic_acid_answers_from_its_whole_table(self):
        self.assertAlmostEqual(d.ph("H3PO4", 0.1), 1.622, places=3)

    def test_an_unknown_name_is_refused(self):
        with self.assertRaises(KeyError):
            d.ph("unobtainium", 0.1)


if __name__ == "__main__":
    unittest.main()
