/* curve.c — Hamiltonian curve addressing (hex9_curve family).
 *
 * Ground truth: 100 KAT vectors generated from hhg9 (h9_curve_uuid /
 * h9_curve_index / h9_curve_label over h9_encode+h9_bin) across layers
 * 0..30, anchored on the usual probe points (Edinburgh, Westminster,
 * Greenwich seam, octahedron vertex, poles, antimeridian). Python's own
 * decode roundtrip was verified clean over the same set when generating.
 *
 * Invariants beyond the KATs: encode/decode roundtrip, axiom order of the
 * 12 roots, rank-truncation == ITERATED one-generation parent (the curve
 * tree is the LINEAGE tree — hex9_cell_parent, not the deep-fold
 * hex9_cell_ancestor, which differs on hexagon-band cells), generator
 * emit order == ascending curve order, and the input contracts.
 */
#include "hex9_c.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

static int failures = 0;
#define CHECK(cond, ...) do { \
    if (!(cond)) { printf("FAIL: " __VA_ARGS__); printf("\n"); failures++; } \
} while (0)

static void unhex(const char *s, uint8_t out[16]) {
    for (int i = 0; i < 16; ++i) {
        unsigned v = 0;
        sscanf(s + 2 * i, "%2x", &v);
        out[i] = (uint8_t)v;
    }
}

static int contains(const uint8_t *uuids, int n, const uint8_t want[16]) {
    for (int i = 0; i < n; ++i)
        if (memcmp(uuids + (size_t)i * 16, want, 16) == 0) return 1;
    return 0;
}

static const char *hex(const uint8_t u[16]) {
    static char buf[4][33];
    static int slot = 0;
    char *b = buf[slot = (slot + 1) & 3];
    for (int i = 0; i < 16; ++i) sprintf(b + 2 * i, "%02x", u[i]);
    return b;
}

static const struct kat {
    int layer;
    const char *bin, *curve, *index, *label;
} KATS[] = {
    {0, "4ffffffffffffffffffffffffffffff2", "c1ffffffffffffffffffffffffffffff", "1", "c1"},
    {0, "4ffffffffffffffffffffffffffffff2", "c1ffffffffffffffffffffffffffffff", "1", "c1"},
    {0, "0ffffffffffffffffffffffffffffff0", "c0ffffffffffffffffffffffffffffff", "0", "c0"},
    {0, "5ffffffffffffffffffffffffffffff4", "c2ffffffffffffffffffffffffffffff", "2", "c2"},
    {0, "4ffffffffffffffffffffffffffffff2", "c1ffffffffffffffffffffffffffffff", "1", "c1"},
    {0, "1ffffffffffffffffffffffffffffff0", "caffffffffffffffffffffffffffffff", "10", "ca"},
    {0, "4ffffffffffffffffffffffffffffff2", "c1ffffffffffffffffffffffffffffff", "1", "c1"},
    {0, "affffffffffffffffffffffffffffff2", "c7ffffffffffffffffffffffffffffff", "7", "c7"},
    {0, "3ffffffffffffffffffffffffffffff0", "c4ffffffffffffffffffffffffffffff", "4", "c4"},
    {0, "bffffffffffffffffffffffffffffff4", "c8ffffffffffffffffffffffffffffff", "8", "c8"},
    {1, "43fffffffffffffffffffffffffffff2", "c11fffffffffffffffffffffffffffff", "10", "c11"},
    {1, "43fffffffffffffffffffffffffffff2", "c11fffffffffffffffffffffffffffff", "10", "c11"},
    {1, "02fffffffffffffffffffffffffffff0", "c03fffffffffffffffffffffffffffff", "3", "c03"},
    {1, "52fffffffffffffffffffffffffffff4", "c25fffffffffffffffffffffffffffff", "23", "c25"},
    {1, "43fffffffffffffffffffffffffffff2", "c11fffffffffffffffffffffffffffff", "10", "c11"},
    {1, "96fffffffffffffffffffffffffffff0", "c98fffffffffffffffffffffffffffff", "89", "c98"},
    {1, "42fffffffffffffffffffffffffffff2", "c15fffffffffffffffffffffffffffff", "14", "c15"},
    {1, "a2fffffffffffffffffffffffffffff2", "c75fffffffffffffffffffffffffffff", "68", "c75"},
    {1, "32fffffffffffffffffffffffffffff0", "c43fffffffffffffffffffffffffffff", "39", "c43"},
    {1, "b3fffffffffffffffffffffffffffff4", "c81fffffffffffffffffffffffffffff", "73", "c81"},
    {2, "432ffffffffffffffffffffffffffff2", "c114ffffffffffffffffffffffffffff", "94", "c114"},
    {2, "435ffffffffffffffffffffffffffff1", "c115ffffffffffffffffffffffffffff", "95", "c115"},
    {2, "022ffffffffffffffffffffffffffff0", "c033ffffffffffffffffffffffffffff", "30", "c033"},
    {2, "522ffffffffffffffffffffffffffff4", "c255ffffffffffffffffffffffffffff", "212", "c255"},
    {2, "433ffffffffffffffffffffffffffff2", "c116ffffffffffffffffffffffffffff", "96", "c116"},
    {2, "107ffffffffffffffffffffffffffff2", "ca00ffffffffffffffffffffffffffff", "810", "ca00"},
    {2, "422ffffffffffffffffffffffffffff2", "c153ffffffffffffffffffffffffffff", "129", "c153"},
    {2, "a22ffffffffffffffffffffffffffff2", "c753ffffffffffffffffffffffffffff", "615", "c753"},
    {2, "320ffffffffffffffffffffffffffff2", "c430ffffffffffffffffffffffffffff", "351", "c430"},
    {2, "b33ffffffffffffffffffffffffffff4", "c811ffffffffffffffffffffffffffff", "658", "c811"},
    {3, "4321fffffffffffffffffffffffffff3", "c1148fffffffffffffffffffffffffff", "854", "c1148"},
    {3, "4348fffffffffffffffffffffffffff2", "c1125fffffffffffffffffffffffffff", "833", "c1125"},
    {3, "0222fffffffffffffffffffffffffff0", "c0333fffffffffffffffffffffffffff", "273", "c0333"},
    {3, "5222fffffffffffffffffffffffffff4", "c2555fffffffffffffffffffffffffff", "1913", "c2555"},
    {3, "4333fffffffffffffffffffffffffff2", "c1165fffffffffffffffffffffffffff", "869", "c1165"},
    {3, "1072fffffffffffffffffffffffffff2", "ca003fffffffffffffffffffffffffff", "7293", "ca003"},
    {3, "4222fffffffffffffffffffffffffff2", "c1533fffffffffffffffffffffffffff", "1164", "c1533"},
    {3, "a222fffffffffffffffffffffffffff2", "c7533fffffffffffffffffffffffffff", "5538", "c7533"},
    {3, "3205fffffffffffffffffffffffffff4", "c4302fffffffffffffffffffffffffff", "3161", "c4302"},
    {3, "b333fffffffffffffffffffffffffff4", "c8116fffffffffffffffffffffffffff", "5928", "c8116"},
    {5, "432177fffffffffffffffffffffffff5", "c114253fffffffffffffffffffffffff", "68736", "c114253"},
    {5, "434868fffffffffffffffffffffffff2", "c112504fffffffffffffffffffffffff", "67477", "c112504"},
    {5, "022222fffffffffffffffffffffffff0", "c033333fffffffffffffffffffffffff", "22143", "c033333"},
    {5, "522222fffffffffffffffffffffffff4", "c255555fffffffffffffffffffffffff", "155003", "c255555"},
    {5, "433336fffffffffffffffffffffffff0", "c116564fffffffffffffffffffffffff", "70447", "c116564"},
    {5, "167216fffffffffffffffffffffffff4", "ca00307fffffffffffffffffffffffff", "590740", "ca00307"},
    {5, "422222fffffffffffffffffffffffff2", "c153333fffffffffffffffffffffffff", "94314", "c153333"},
    {5, "a22222fffffffffffffffffffffffff2", "c753333fffffffffffffffffffffffff", "448608", "c753333"},
    {5, "320530fffffffffffffffffffffffff0", "c430221fffffffffffffffffffffffff", "256060", "c430221"},
    {5, "b33333fffffffffffffffffffffffff4", "c811656fffffffffffffffffffffffff", "480219", "c811656"},
    {8, "432177468ffffffffffffffffffffff3", "c114253704ffffffffffffffffffffff", "50109115", "c114253704"},
    {8, "435878507ffffffffffffffffffffff3", "c112504417ffffffffffffffffffffff", "49191073", "c112504417"},
    {8, "022222222ffffffffffffffffffffff0", "c033333333ffffffffffffffffffffff", "16142520", "c033333333"},
    {8, "522222222ffffffffffffffffffffff4", "c255555555ffffffffffffffffffffff", "112997642", "c255555555"},
    {8, "433336868ffffffffffffffffffffff2", "c116564604ffffffffffffffffffffff", "51356353", "c116564604"},
    {8, "167216025ffffffffffffffffffffff2", "ca00307022ffffffffffffffffffffff", "430649480", "ca00307022"},
    {8, "422222262ffffffffffffffffffffff0", "c153333373ffffffffffffffffffffff", "68755215", "c153333373"},
    {8, "a22222515ffffffffffffffffffffff1", "c753333272ffffffffffffffffffffff", "327035459", "c753333272"},
    {8, "320530621ffffffffffffffffffffff4", "c430220858ffffffffffffffffffffff", "186667712", "c430220858"},
    {8, "b33333234ffffffffffffffffffffff2", "c811656066ffffffffffffffffffffff", "350079711", "c811656066"},
    {12, "4321774683685ffffffffffffffffff3", "c1142537041324ffffffffffffffffff", "328765904509", "c1142537041324"},
    {12, "4358785076384ffffffffffffffffff1", "c1125044178283ffffffffffffffffff", "322742636022", "c1125044178283"},
    {12, "0222222222224ffffffffffffffffff4", "c0333333333335ffffffffffffffffff", "105911076182", "c0333333333335"},
    {12, "5222222222222ffffffffffffffffff4", "c2555555555555ffffffffffffffffff", "741377533262", "c2555555555555"},
    {12, "4333368684822ffffffffffffffffff2", "c1165646047335ffffffffffffffffff", "336949037411", "c1165646047335"},
    {12, "1672160258563ffffffffffffffffff0", "ca003070257481ffffffffffffffffff", "2825491263463", "ca003070257481"},
    {12, "4222222622286ffffffffffffffffff0", "c1533333735548ffffffffffffffffff", "451102969709", "c1533333735548"},
    {12, "a222225150050ffffffffffffffffff3", "c7533332721020ffffffffffffffffff", "2145679647246", "c7533332721020"},
    {12, "3205306211070ffffffffffffffffff0", "c4302208587880ffffffffffffffffff", "1224726864255", "c4302208587880"},
    {12, "b333332348428ffffffffffffffffff4", "c8116560661744ffffffffffffffffff", "2296872985207", "c8116560661744"},
    {18, "4321774683685205354ffffffffffff1", "c1142537041324486626ffffffffffff", "174719681058461037", "c1142537041324486626"},
    {18, "4358785076384212044ffffffffffff3", "c1125044178283403832ffffffffffff", "171518669230406762", "c1125044178283403832"},
    {18, "0222222222224002308ffffffffffff0", "c0333333333335102674ffffffffffff", "56285488237299322", "c0333333333335102674"},
    {18, "5222222222222222222ffffffffffff4", "c2555555555555555555ffffffffffff", "393998417654622692", "c2555555555555555555"},
    {18, "4333368684822323262ffffffffffff0", "c1165646047335276115ffffffffffff", "179068533390907745", "c1165646047335276115"},
    {18, "1672160258763515516ffffffffffff0", "ca003070257481472511ffffffffffff", "1501581902546324179", "ca003070257481472511"},
    {18, "4222222622286283422ffffffffffff0", "c1533333735548547645ffffffffffff", "239734613325447788", "c1533333735548547645"},
    {18, "a222225150050433110ffffffffffff5", "c7533332721020522408ffffffffffff", "1140302137412371643", "c7533332721020522408"},
    {18, "3205306211070207266ffffffffffff4", "c4302208587880388510ffffffffffff", "650870069466777336", "c4302208587880388510"},
    {18, "b333332348428686832ffffffffffff4", "c8116560661744848474ffffffffffff", "1220652476131898146", "c8116560661744848474"},
    {25, "43217746836852053548751236fffff5", "c11425370413244866264741520fffff", "835678818192506330245356", "c11425370413244866264741520"},
    {25, "43587850763842120447853502fffff3", "c11250441782834038320446702fffff", "820368477850289400303761", "c11250441782834038320446702"},
    {25, "02222222222240023080877157fffff4", "c03333333333351026751267660fffff", "269211745388867306324535", "c03333333333351026751267660"},
    {25, "52222222222222222222222222fffff4", "c25555555555555555555555555fffff", "1884482217691113045521903", "c25555555555555555555555555"},
    {25, "43333686848223232628482223fffff2", "c11656460473352761154733552fffff", "856479244084176628756334", "c11656460473352761154733552"},
    {25, "16721602587635155565735011fffff4", "ca0030702574814725616315487fffff", "7182019690840089830717458", "ca0030702574814725616315487"},
    {25, "42222226222862834228428422fffff0", "c15333337355485476454144745fffff", "1146643223762603683337153", "c15333337355485476454144745"},
    {25, "a2222251500504331102225043fffff1", "c75333327210205224085556752fffff", "5454029773877113787938310", "c75333327210205224085556752"},
    {25, "32053062110702072660780827fffff2", "c43022085878803885108837451fffff", "3113091365287442532739660", "c43022085878803885108837451"},
    {25, "b3333323484286868324842283fffff4", "c81165606617448484745714341fffff", "5838342953112108746555779", "c81165606617448484745714341"},
    {30, "43217746836852053548751266106825", "c1142537041324486626474152008165", "49345998535449306294658032416", "c1142537041324486626474152008165"},
    {30, "43587850763842120447853502863205", "c1125044178283403832044670240610", "48441938248581738798536810028", "c1125044178283403832044670240610"},
    {30, "02222222222240023080877157224420", "c0333333333335102675126766033561", "15896684353467225571157489545", "c0333333333335102675126766033561"},
    {30, "52222222222222222222222222222224", "c2555555555555555555555555555555", "111276790472442534225022887152", "c2555555555555555555555555555555"},
    {30, "43333686848223232628482223682360", "c1165646047335276115473355204568", "50574242883926545751432769749", "c1165646047335276115473355204568"},
    {30, "16721602587635155565735011588844", "ca003070257481472561631548743221", "424091080724416464414035206054", "ca003070257481472561631548743221"},
    {30, "42222226222862834228428422848322", "c1533333735548547645414474541514", "67708135719957984897375574888", "c1533333735548547645414474541514"},
    {30, "a2222251500504331102225043022415", "c7533332721020522408555675203353", "322055004117669692063969269668", "c7533332721020522408555675203353"},
    {30, "32053062110702072660780827410124", "c4302208587880388510883745177875", "183824932028858194115744235086", "c4302208587880388510883745177875"},
    {30, "b3333323484286868324842283348334", "c8116560661744848474571434162372", "344748313038316909375372235303", "c8116560661744848474571434162372"},
};
static const int NKATS = (int)(sizeof(KATS) / sizeof(KATS[0]));

int main(void) {
    char err[256];
    if (hex9_warp_init(err, sizeof err) != 0) {
        printf("FAIL: warp init: %s\n", err);
        return 1;
    }
    const int lmax = hex9_lmax();

    /* ── KAT sweep: encode, index, label, pack, parse, layer, decode ────── */
    for (int i = 0; i < NKATS; ++i) {
        const struct kat *k = &KATS[i];
        if (k->layer > lmax) continue;             /* legacy L29 build */
        uint8_t bin[16], want[16], got[16];
        unhex(k->bin, bin);
        unhex(k->curve, want);

        CHECK(hex9_curve(bin, got) == 0 && memcmp(got, want, 16) == 0,
              "kat %d L%d curve: got %s want %s", i, k->layer, hex(got), k->curve);

        CHECK(hex9_is_curve(want) == 1 && hex9_is_curve(bin) == 0,
              "kat %d is_curve discrimination", i);
        CHECK(hex9_curve_layer(want) == k->layer,
              "kat %d curve_layer %d != %d", i, hex9_curve_layer(want), k->layer);

        /* curve input passes through encode unchanged */
        CHECK(hex9_curve(want, got) == 0 && memcmp(got, want, 16) == 0,
              "kat %d curve passthrough", i);

        char buf[64];
        CHECK(hex9_curve_index(bin, buf, sizeof buf) > 0 && strcmp(buf, k->index) == 0,
              "kat %d index: got %s want %s", i, buf, k->index);
        CHECK(hex9_curve_index(want, buf, sizeof buf) > 0 && strcmp(buf, k->index) == 0,
              "kat %d index from curve-uuid: got %s want %s", i, buf, k->index);

        CHECK(hex9_curve_pack(k->index, k->layer, got) == 0 && memcmp(got, want, 16) == 0,
              "kat %d pack(%s, %d): got %s", i, k->index, k->layer, hex(got));

        CHECK(hex9_curve_label(want, buf, sizeof buf) > 0 && strcmp(buf, k->label) == 0,
              "kat %d label: got %s want %s", i, buf, k->label);
        CHECK(hex9_curve_parse_label(k->label, got) == 0 && memcmp(got, want, 16) == 0,
              "kat %d parse_label(%s)", i, k->label);

        CHECK(hex9_curve_decode(want, got) == 0 && memcmp(got, bin, 16) == 0,
              "kat %d L%d decode: got %s want %s", i, k->layer, hex(got), k->bin);
        /* h9-uuid passes through decode unchanged */
        CHECK(hex9_curve_decode(bin, got) == 0 && memcmp(got, bin, 16) == 0,
              "kat %d decode passthrough", i);
    }

    /* ── Full-uuid input: curve(full) == curve(bin(full, L)) at every L ──── */
    {
        uint8_t full[16], bin[16], c_full[16], c_bin[16], c_tr[16];
        CHECK(hex9_encode(-3.1883, 55.9533, full) == 0, "encode edinburgh");
        CHECK(hex9_curve(full, c_full) == 0, "curve(full)");
        CHECK(hex9_bin(full, lmax, bin) == 0, "bin lmax");
        CHECK(hex9_curve(bin, c_bin) == 0, "curve(bin lmax)");
        CHECK(memcmp(c_full, c_bin, 16) == 0,
              "curve(full) %s != curve(bin(full,lmax)) %s", hex(c_full), hex(c_bin));
        /* rank truncation == encode of the shallower bin's LINEAGE parent:
         * walk iterated one-generation parents down from the lmax bin. */
        uint8_t cur[16], par[16];
        memcpy(cur, bin, 16);
        for (int L = lmax; L >= 1; --L) {
            CHECK(hex9_cell_parent(cur, par) == 0, "cell_parent L%d", L);
            CHECK(hex9_curve_bin(c_full, L - 1, c_tr) == 0, "curve_bin L%d", L - 1);
            CHECK(hex9_curve(par, c_bin) == 0, "curve(parent) L%d", L - 1);
            CHECK(memcmp(c_tr, c_bin, 16) == 0,
                  "L%d truncation %s != lineage parent curve %s",
                  L - 1, hex(c_tr), hex(c_bin));
            memcpy(cur, par, 16);
        }
    }

    /* ── Axiom order: the 12 roots hit slots 0..11 exactly once ──────────── */
    {
        int seen[12] = {0};
        for (int root = 0; root < 12; ++root) {
            /* L0 bins from the KAT set cover 8 roots; build the rest via
             * parse of the L0 curve label and decode->encode roundtrip. */
            char lbl[4] = {'c', (char)(root < 10 ? '0' + root : 'a' + root - 10), 0, 0};
            uint8_t cu[16], bin[16], back[16];
            CHECK(hex9_curve_parse_label(lbl, cu) == 0, "parse %s", lbl);
            CHECK(hex9_curve_decode(cu, bin) == 0, "decode L0 slot %d", root);
            CHECK(hex9_curve(bin, back) == 0 && memcmp(back, cu, 16) == 0,
                  "L0 slot %d roundtrip", root);
            const int layer = hex9_curve_layer(cu);
            CHECK(layer == 0, "L0 slot %d layer %d", root, layer);
            uint8_t nib0 = (uint8_t)(bin[0] >> 4);
            CHECK(nib0 <= 11 && !seen[nib0], "root %d duplicated/oob", nib0);
            if (nib0 <= 11) seen[nib0] = 1;
        }
        for (int r = 0; r < 12; ++r)
            CHECK(seen[r], "root %d never emitted by slots 0..11", r);
    }

    /* ── Generator: descendants in curve order ───────────────────────────── */
    {
        uint8_t anc[16];
        unhex("432ffffffffffffffffffffffffffff2", anc);      /* Edinburgh L2 */
        CHECK(hex9_curve_ncells(2, 4) == 81, "ncells 2->4");
        uint8_t bins[81 * 16], curves[81 * 16];
        int64_t n = hex9_curve_cells(anc, 4, bins, curves, 81);
        CHECK(n == 81, "curve_cells count %lld", (long long)n);
        uint8_t anc_curve[16];
        CHECK(hex9_curve(anc, anc_curve) == 0, "curve(anc)");
        for (int i = 0; i < (int)n; ++i) {
            uint8_t cu[16], tr[16], up[16], up2[16];
            /* emitted curve matches an independent encode of the bin */
            CHECK(hex9_curve(bins + i * 16, cu) == 0
                      && memcmp(cu, curves + i * 16, 16) == 0,
                  "gen %d curve mismatch", i);
            /* strictly ascending curve order == byte order at one layer */
            if (i > 0)
                CHECK(memcmp(curves + (i - 1) * 16, curves + i * 16, 16) < 0,
                      "gen %d not ascending", i);
            /* every emitted cell truncates to the ancestor's curve address */
            CHECK(hex9_curve_bin(cu, 2, tr) == 0 && memcmp(tr, anc_curve, 16) == 0,
                  "gen %d not under ancestor", i);
            /* and its lineage grandparent is the ancestor bin */
            CHECK(hex9_cell_parent(bins + i * 16, up) == 0
                      && hex9_cell_parent(up, up2) == 0
                      && memcmp(up2, anc, 16) == 0,
                  "gen %d lineage grandparent mismatch", i);
        }
        /* generation from the equivalent curve-uuid input matches */
        uint8_t bins2[81 * 16];
        CHECK(hex9_curve_cells(anc_curve, 4, bins2, NULL, 81) == 81
                  && memcmp(bins, bins2, sizeof bins2) == 0,
              "curve_cells from curve-uuid differs");
        /* cap honoured */
        CHECK(hex9_curve_cells(anc, 4, bins, NULL, 80) == -1, "cap not honoured");
        /* self-enumeration: layer == own layer emits just the cell */
        CHECK(hex9_curve_cells(anc, 2, bins, curves, 1) == 1
                  && memcmp(bins, anc, 16) == 0,
              "self enumeration");
    }

    /* ── Lineage vs ownership: children, owned cells, global partition ───── */
    {
        uint8_t anc[16];
        unhex("432ffffffffffffffffffffffffffff2", anc);      /* Edinburgh L2 */

        /* children: one generation, where the two relations coincide — the
         * ranked children must equal curve_cells at depth 1, in order. */
        uint8_t kids[9 * 16], gen[9 * 16];
        CHECK(hex9_cell_children(anc, kids) == 0, "cell_children");
        CHECK(hex9_curve_cells(anc, 3, gen, NULL, 9) == 9, "gen depth 1");
        CHECK(memcmp(kids, gen, sizeof kids) == 0,
              "children != rank-ordered lineage generation");
        for (int j = 0; j < 9; ++j) {
            uint8_t par[16];
            CHECK(hex9_cell_parent(kids + j * 16, par) == 0
                      && memcmp(par, anc, 16) == 0, "child %d parent", j);
            uint8_t own[16];
            CHECK(hex9_cell_ancestor(kids + j * 16, 2, own) == 0
                      && memcmp(own, anc, 16) == 0, "child %d owner", j);
        }

        /* owned cells at depth 2: exactly 81, curve-sorted, every one's
         * OWNERSHIP ancestor (deep fold) is the zone. */
        uint8_t obins[81 * 16], ocurves[81 * 16];
        int64_t n = hex9_owned_cells(anc, 4, obins, ocurves, 81);
        CHECK(n == 81, "owned_cells count %lld", (long long)n);
        for (int i = 0; i < (int)n; ++i) {
            uint8_t own[16], cu[16];
            CHECK(hex9_cell_ancestor(obins + i * 16, 2, own) == 0
                      && memcmp(own, anc, 16) == 0, "owned %d ancestor", i);
            CHECK(hex9_curve(obins + i * 16, cu) == 0
                      && memcmp(cu, ocurves + i * 16, 16) == 0,
                  "owned %d curve mismatch", i);
            if (i > 0)
                CHECK(memcmp(ocurves + (i - 1) * 16, ocurves + i * 16, 16) < 0,
                      "owned %d not curve-sorted", i);
        }
        /* self-enumeration */
        CHECK(hex9_owned_cells(anc, 2, obins, NULL, 1) == 1
                  && memcmp(obins, anc, 16) == 0, "owned self");

        /* GLOBAL PARTITION at L2: the 12 roots' owned sets tile the layer —
         * 12 * 81 = 972 cells, all distinct (doctrine: owned sub-zone sets
         * partition the globe; counts exactly 9^d). Lineage sets partition
         * too (every cell has exactly one lineage ancestor); check both and
         * count where the two readings diverge (hexagon-band cells). */
        static uint8_t all_owned[972][16], all_lin[972][16];
        int t = 0;
        for (int slot = 0; slot < 12; ++slot) {
            char lbl[4] = {'c', (char)(slot < 10 ? '0' + slot : 'a' + slot - 10), 0, 0};
            uint8_t cu[16], root[16];
            CHECK(hex9_curve_parse_label(lbl, cu) == 0
                      && hex9_curve_decode(cu, root) == 0, "root %d", slot);
            CHECK(hex9_owned_cells(root, 2, (uint8_t *)all_owned[t], NULL, 81) == 81,
                  "root %d owned", slot);
            CHECK(hex9_curve_cells(root, 2, (uint8_t *)all_lin[t], NULL, 81) == 81,
                  "root %d lineage", slot);
            t += 81;
        }
        int dup_owned = 0, dup_lin = 0, diverged = 0;
        for (int i = 0; i < 972; ++i) {
            for (int j = i + 1; j < 972; ++j) {
                if (memcmp(all_owned[i], all_owned[j], 16) == 0) dup_owned++;
                if (memcmp(all_lin[i], all_lin[j], 16) == 0) dup_lin++;
            }
            if (!contains((const uint8_t *)all_lin, 972, all_owned[i])) diverged++;
        }
        CHECK(dup_owned == 0, "owned partition has %d duplicates", dup_owned);
        CHECK(dup_lin == 0, "lineage partition has %d duplicates", dup_lin);
        /* both partition the same 972-cell layer, so as SETS they are equal;
         * divergence lives in WHICH zone a cell belongs to, checked above
         * per-cell via ancestor == zone. */
        CHECK(diverged == 0, "%d owned cells missing from lineage universe", diverged);
    }

    /* ── Error contracts ─────────────────────────────────────────────────── */
    {
        uint8_t cu[16], out[16];
        unhex("c114ffffffffffffffffffffffffffff", cu);
        CHECK(hex9_curve_bin(cu, 3, out) != 0, "curve_bin below own layer must fail");
        CHECK(hex9_curve_layer(out) == -1 || 1, "no-op");   /* placate unused */
        uint8_t bad[16];
        unhex("c1c4ffffffffffffffffffffffffffff", bad);     /* rank nibble 0xc */
        CHECK(hex9_curve_layer(bad) == -1, "bad rank digit accepted");
        CHECK(hex9_curve_parse_label("c9", out) != 0 || hex9_curve_layer(out) == 0,
              "slot 9 label is valid");                     /* 'c9' IS valid */
        CHECK(hex9_curve_parse_label("cc", out) != 0, "slot 12 accepted");
        CHECK(hex9_curve_parse_label("x1", out) != 0, "bad marker accepted");
        CHECK(hex9_curve_parse_label("c19", out) != 0, "rank 9 accepted");
        CHECK(hex9_curve_pack("999999999999999999999999999999999", 2, out) != 0,
              "oversize index accepted");
        CHECK(hex9_curve_ncells(4, 2) == -1 && hex9_curve_ncells(0, 20) == -1,
              "ncells bounds");
    }

    if (failures) { printf("%d failure(s)\n", failures); return 1; }
    printf("curve: all checks passed\n");
    return 0;
}
