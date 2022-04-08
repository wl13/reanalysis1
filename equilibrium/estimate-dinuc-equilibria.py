### IMPORTS ###

import os
import numpy as np
import csv
import random
import scipy.stats

### SOURCE FILE ###

source = 'path/to/mutational-matrix.tsv'
out_source = 'path/to/output.tsv'

### MAIN ###

def main():

    #Extract mutational data from the input file
    raw_file = open(source).read()
    split = raw_file.split('\n')
    counts = [float(i.split('\t')[2]) for i in split[1:] if i != '']
    totals = [float(i.split('\t')[3]) for i in split[1:] if i != '']

    #Calculate mutational frequencies
    freqs = [x/y for x, y in zip(counts, totals)]

    #Use mutational frequencies to estimate mutational equilibrium
    aa_e, ac_e, ag_e, at_e, ca_e, cc_e, cg_e, ct_e, ga_e, gc_e, gg_e, gt_e, ta_e, tc_e, tg_e, tt_e, predicted_gc_e = get_dinuc_eq(freqs)

    #Bootstrap

    #Obtain a list of all mutations
    ref = [i.split('\t')[0] for i in split[1:] if i != '']
    alt = [i.split('\t')[1] for i in split[1:] if i != '']
    mutations = [str(x) + str(y) for x, y in zip(ref, alt)]
    sampling_list = []
    for i,n in enumerate(mutations):
        name = n
        count = counts[i]
        no = 0
        while no < count:
            sampling_list.append(name)
            no += 1

    #Create empty lists to contain bootstrap estimates
    aa = []
    ac = []
    ag = []
    at = []
    ca = []
    cc = []
    cg = []
    ct = []
    ga = []
    gc = []
    gg = []
    gt = []
    ta = []
    tc = []
    tg = []
    tt = []
    gc_content = []

    #Resample with replacement to receive a bootstrap, then estimate equilibrium for each dinucleotide and add to above lists
    counter = 0
    while counter < 100: #edit to determine the number of bootstraps to be completed
        sample = random.choices(sampling_list, k=len(sampling_list))
        new_counts = []
        for mutation in mutations:
            count = sample.count(mutation)
            new_counts.append(count)
        frequencies = [x/y for x, y in zip(new_counts, totals)]
        aa_eq, ac_eq, ag_eq, at_eq, ca_eq, cc_eq, cg_eq, ct_eq, ga_eq, gc_eq, gg_eq, gt_eq, ta_eq, tc_eq, tg_eq, tt_eq, predicted_gc = get_dinuc_eq(frequencies)
        aa.append(aa_eq)
        ac.append(ac_eq)
        ag.append(ag_eq)
        at.append(at_eq)
        ca.append(ca_eq)
        cc.append(cc_eq)
        cg.append(cg_eq)
        ct.append(ct_eq)
        ga.append(ga_eq)
        gc.append(gc_eq)
        gg.append(gg_eq)
        gt.append(gt_eq)
        ta.append(ta_eq)
        tc.append(tc_eq)
        tg.append(tg_eq)
        tt.append(tt_eq)
        gc_content.append(predicted_gc)
        counter += 1

    #Build output file
    csv_total = []
    headers = ['dinucleotide', 'estimate', 'bootstrapped_mean', 'bootstrapped_lower_95_ci', 'bootstrapped_upper_95_ci']

    mean_aa, lower_aa, upper_aa = mean_confidence_interval(aa, confidence=0.95)
    print ('Estimate, mean, lower, upper:')
    print ('AA:', aa_e, mean_aa, lower_aa, upper_aa)
    csv_total.append(['AA', aa_e, mean_aa, lower_aa, upper_aa])

    mean_ac, lower_ac, upper_ac = mean_confidence_interval(ac, confidence=0.95)
    print ('AC:', ac_e, mean_ac, lower_ac, upper_ac)
    csv_total.append(['AC', ac_e, mean_ac, lower_ac, upper_ac])

    mean_ag, lower_ag, upper_ag = mean_confidence_interval(ag, confidence=0.95)
    print ('AG:', ag_e, mean_ag, lower_ag, upper_ag)
    csv_total.append(['AG', ag_e, mean_ag, lower_ag, upper_ag])

    mean_at, lower_at, upper_at = mean_confidence_interval(at, confidence=0.95)
    print ('AT:', at_e, mean_at, lower_at, upper_at)
    csv_total.append(['AT', at_e, mean_at, lower_at, upper_at])

    mean_ca, lower_ca, upper_ca = mean_confidence_interval(ca, confidence=0.95)
    print ('CA:', ca_e, mean_ca, lower_ca, upper_ca)
    csv_total.append(['CA', ca_e, mean_ca, lower_ca, upper_ca])

    mean_cc, lower_cc, upper_cc = mean_confidence_interval(cc, confidence=0.95)
    print ('CC:', cc_e, mean_cc, lower_cc, upper_cc)
    csv_total.append(['CC', cc_e, mean_cc, lower_cc, upper_cc])

    mean_cg, lower_cg, upper_cg = mean_confidence_interval(cg, confidence=0.95)
    print ('CG:', cg_e, mean_cg, lower_cg, upper_cg)
    csv_total.append(['CG', cg_e, mean_cg, lower_cg, upper_cg])

    mean_ct, lower_ct, upper_ct = mean_confidence_interval(ct, confidence=0.95)
    print ('CT:', ct_e, mean_ct, lower_ct, upper_ct)
    csv_total.append(['CT', ct_e, mean_ct, lower_ct, upper_ct])

    mean_ga, lower_ga, upper_ga = mean_confidence_interval(ga, confidence=0.95)
    print ('GA:', ga_e, mean_ga, lower_ga, upper_ga)
    csv_total.append(['GA', ga_e, mean_ga, lower_ga, upper_ga])

    mean_gc, lower_gc, upper_gc = mean_confidence_interval(gc, confidence=0.95)
    print ('GC:', gc_e, mean_gc, lower_gc, upper_gc)
    csv_total.append(['GC', gc_e, mean_gc, lower_gc, upper_gc])

    mean_gg, lower_gg, upper_gg = mean_confidence_interval(gg, confidence=0.95)
    print ('GG:', gg_e, mean_gg, lower_gg, upper_gg)
    csv_total.append(['GG', gg_e, mean_gg, lower_gg, upper_gg])

    mean_gt, lower_gt, upper_gt = mean_confidence_interval(gt, confidence=0.95)
    print ('GT:', gt_e, mean_gt, lower_gt, upper_gt)
    csv_total.append(['GT', gt_e, mean_gt, lower_gt, upper_gt])

    mean_ta, lower_ta, upper_ta = mean_confidence_interval(ta, confidence=0.95)
    print ('TA:', ta_e, mean_ta, lower_ta, upper_ta)
    csv_total.append(['TA', ta_e, mean_ta, lower_ta, upper_ta])

    mean_tc, lower_tc, upper_tc = mean_confidence_interval(tc, confidence=0.95)
    print ('TC:', tc_e, mean_tc, lower_tc, upper_tc)
    csv_total.append(['TC', tc_e, mean_tc, lower_tc, upper_tc])

    mean_tg, lower_tg, upper_tg = mean_confidence_interval(tg, confidence=0.95)
    print ('TG:', tg_e, mean_tg, lower_tg, upper_tg)
    csv_total.append(['TG', tg_e, mean_tg, lower_tg, upper_tg])

    mean_tt, lower_tt, upper_tt = mean_confidence_interval(tt, confidence=0.95)
    print ('TT:', tt_e, mean_tt, lower_tt, upper_tt)
    csv_total.append(['TT', tt_e, mean_tt, lower_tt, upper_tt])

    mean_gc_c, lower_gc_c, upper_gc_c = mean_confidence_interval(gc_content, confidence=0.95)
    print ('GC content:', predicted_gc_e, mean_gc_c, lower_gc_c, upper_gc_c)

    #Write output to predefined source
    filepath = os.path.abspath(out_source)
    with open(filepath, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(i for i in headers)
        for j in csv_total:
            writer.writerow(j)

### FUNCTIONS ###

def get_dinuc_eq(data):

    #First extract data for ease later

    #AA
    aa_ac = float(data[0])
    aa_ag = float(data[1])
    aa_at = float(data[2])
    aa_ca = float(data[3])
    aa_cc = float(data[4])
    aa_cg = float(data[5])
    aa_ct = float(data[6])
    aa_ga = float(data[7])
    aa_gc = float(data[8])
    aa_gg = float(data[9])
    aa_gt = float(data[10])
    aa_ta = float(data[11])
    aa_tc = float(data[12])
    aa_tg = float(data[13])
    aa_tt = float(data[14])

    #AC
    ac_aa = float(data[15])
    ac_ag = float(data[16])
    ac_at = float(data[17])
    ac_ca = float(data[18])
    ac_cc = float(data[19])
    ac_cg = float(data[20])
    ac_ct = float(data[21])
    ac_ga = float(data[22])
    ac_gc = float(data[23])
    ac_gg = float(data[24])
    ac_gt = float(data[25])
    ac_ta = float(data[26])
    ac_tc = float(data[27])
    ac_tg = float(data[28])
    ac_tt = float(data[29])

    #AG
    ag_aa = float(data[30])
    ag_ac = float(data[31])
    ag_at = float(data[32])
    ag_ca = float(data[33])
    ag_cc = float(data[34])
    ag_cg = float(data[35])
    ag_ct = float(data[36])
    ag_ga = float(data[37])
    ag_gc = float(data[38])
    ag_gg = float(data[39])
    ag_gt = float(data[40])
    ag_ta = float(data[41])
    ag_tc = float(data[42])
    ag_tg = float(data[43])
    ag_tt = float(data[44])

    #AT
    at_aa = float(data[45])
    at_ac = float(data[46])
    at_ag = float(data[47])
    at_ca = float(data[48])
    at_cc = float(data[49])
    at_cg = float(data[50])
    at_ct = float(data[51])
    at_ga = float(data[52])
    at_gc = float(data[53])
    at_gg = float(data[54])
    at_gt = float(data[55])
    at_ta = float(data[56])
    at_tc = float(data[57])
    at_tg = float(data[58])
    at_tt = float(data[59])

    #CA
    ca_aa = float(data[60])
    ca_ac = float(data[61])
    ca_ag = float(data[62])
    ca_at = float(data[63])
    ca_cc = float(data[64])
    ca_cg = float(data[65])
    ca_ct = float(data[66])
    ca_ga = float(data[67])
    ca_gc = float(data[68])
    ca_gg = float(data[69])
    ca_gt = float(data[70])
    ca_ta = float(data[71])
    ca_tc = float(data[72])
    ca_tg = float(data[73])
    ca_tt = float(data[74])

    #CC
    cc_aa = float(data[75])
    cc_ac = float(data[76])
    cc_ag = float(data[77])
    cc_at = float(data[78])
    cc_ca = float(data[79])
    cc_cg = float(data[80])
    cc_ct = float(data[81])
    cc_ga = float(data[82])
    cc_gc = float(data[83])
    cc_gg = float(data[84])
    cc_gt = float(data[85])
    cc_ta = float(data[86])
    cc_tc = float(data[87])
    cc_tg = float(data[88])
    cc_tt = float(data[89])

    #CG
    cg_aa = float(data[90])
    cg_ac = float(data[91])
    cg_ag = float(data[92])
    cg_at = float(data[93])
    cg_ca = float(data[94])
    cg_cc = float(data[95])
    cg_ct = float(data[96])
    cg_ga = float(data[97])
    cg_gc = float(data[98])
    cg_gg = float(data[99])
    cg_gt = float(data[100])
    cg_ta = float(data[101])
    cg_tc = float(data[102])
    cg_tg = float(data[103])
    cg_tt = float(data[104])

    #CT
    ct_aa = float(data[105])
    ct_ac = float(data[106])
    ct_ag = float(data[107])
    ct_at = float(data[108])
    ct_ca = float(data[109])
    ct_cc = float(data[110])
    ct_cg = float(data[111])
    ct_ga = float(data[112])
    ct_gc = float(data[113])
    ct_gg = float(data[114])
    ct_gt = float(data[115])
    ct_ta = float(data[116])
    ct_tc = float(data[117])
    ct_tg = float(data[118])
    ct_tt = float(data[119])

    #GA
    ga_aa = float(data[120])
    ga_ac = float(data[121])
    ga_ag = float(data[122])
    ga_at = float(data[123])
    ga_ca = float(data[124])
    ga_cc = float(data[125])
    ga_cg = float(data[126])
    ga_ct = float(data[127])
    ga_gc = float(data[128])
    ga_gg = float(data[129])
    ga_gt = float(data[130])
    ga_ta = float(data[131])
    ga_tc = float(data[132])
    ga_tg = float(data[133])
    ga_tt = float(data[134])

    #GC
    gc_aa = float(data[135])
    gc_ac = float(data[136])
    gc_ag = float(data[137])
    gc_at = float(data[138])
    gc_ca = float(data[139])
    gc_cc = float(data[140])
    gc_cg = float(data[141])
    gc_ct = float(data[142])
    gc_ga = float(data[143])
    gc_gg = float(data[144])
    gc_gt = float(data[145])
    gc_ta = float(data[146])
    gc_tc = float(data[147])
    gc_tg = float(data[148])
    gc_tt = float(data[149])

    #GG
    gg_aa = float(data[150])
    gg_ac = float(data[151])
    gg_ag = float(data[152])
    gg_at = float(data[153])
    gg_ca = float(data[154])
    gg_cc = float(data[155])
    gg_cg = float(data[156])
    gg_ct = float(data[157])
    gg_ga = float(data[158])
    gg_gc = float(data[159])
    gg_gt = float(data[160])
    gg_ta = float(data[161])
    gg_tc = float(data[162])
    gg_tg = float(data[163])
    gg_tt = float(data[164])

    #GT
    gt_aa = float(data[165])
    gt_ac = float(data[166])
    gt_ag = float(data[167])
    gt_at = float(data[168])
    gt_ca = float(data[169])
    gt_cc = float(data[170])
    gt_cg = float(data[171])
    gt_ct = float(data[172])
    gt_ga = float(data[173])
    gt_gc = float(data[174])
    gt_gg = float(data[175])
    gt_ta = float(data[176])
    gt_tc = float(data[177])
    gt_tg = float(data[178])
    gt_tt = float(data[179])

    #TA
    ta_aa = float(data[180])
    ta_ac = float(data[181])
    ta_ag = float(data[182])
    ta_at = float(data[183])
    ta_ca = float(data[184])
    ta_cc = float(data[185])
    ta_cg = float(data[186])
    ta_ct = float(data[187])
    ta_ga = float(data[188])
    ta_gc = float(data[189])
    ta_gg = float(data[190])
    ta_gt = float(data[191])
    ta_tc = float(data[192])
    ta_tg = float(data[193])
    ta_tt = float(data[194])

    #TC
    tc_aa = float(data[195])
    tc_ac = float(data[196])
    tc_ag = float(data[197])
    tc_at = float(data[198])
    tc_ca = float(data[199])
    tc_cc = float(data[200])
    tc_cg = float(data[201])
    tc_ct = float(data[202])
    tc_ga = float(data[203])
    tc_gc = float(data[204])
    tc_gg = float(data[205])
    tc_gt = float(data[206])
    tc_ta = float(data[207])
    tc_tg = float(data[208])
    tc_tt = float(data[209])

    #TG
    tg_aa = float(data[210])
    tg_ac = float(data[211])
    tg_ag = float(data[212])
    tg_at = float(data[213])
    tg_ca = float(data[214])
    tg_cc = float(data[215])
    tg_cg = float(data[216])
    tg_ct = float(data[217])
    tg_ga = float(data[218])
    tg_gc = float(data[219])
    tg_gg = float(data[220])
    tg_gt = float(data[221])
    tg_ta = float(data[222])
    tg_tc = float(data[223])
    tg_tt = float(data[224])

    #TT
    tt_aa = float(data[225])
    tt_ac = float(data[226])
    tt_ag = float(data[227])
    tt_at = float(data[228])
    tt_ca = float(data[229])
    tt_cc = float(data[230])
    tt_cg = float(data[231])
    tt_ct = float(data[232])
    tt_ga = float(data[233])
    tt_gc = float(data[234])
    tt_gg = float(data[235])
    tt_gt = float(data[236])
    tt_ta = float(data[237])
    tt_tc = float(data[238])
    tt_tg = float(data[239])

    #Equation 1
    aa1 = 1 - (1 - aa_ac - aa_ag - aa_at - aa_ca - aa_cc - aa_cg - aa_ct - aa_ga - aa_gc - aa_gg - aa_gt - aa_ta - aa_tc - aa_tg - aa_tt) + tt_aa
    ac1 = tt_aa - ac_aa
    ag1 = tt_aa - ag_aa
    at1 = tt_aa - at_aa
    ca1 = tt_aa - ca_aa
    cc1 = tt_aa - cc_aa
    cg1 = tt_aa - cg_aa
    ct1 = tt_aa - ct_aa
    ga1 = tt_aa - ga_aa
    gc1 = tt_aa - gc_aa
    gg1 = tt_aa - gg_aa
    gt1 = tt_aa - gt_aa
    ta1 = tt_aa - ta_aa
    tc1 = tt_aa - tc_aa
    tg1 = tt_aa - tg_aa
    constant1 = tt_aa
    equation1 = [aa1, ac1, ag1, at1, ca1, cc1, cg1, ct1, ga1, gc1, gg1, gt1, ta1, tc1, tg1]

    #Equation 2
    aa2 = tt_ac - aa_ac
    ac2 = 1 - (1 - ac_aa - ac_ag - ac_at - ac_ca - ac_cc - ac_cg - ac_ct - ac_ga - ac_gc - ac_gg - ac_gt - ac_ta - ac_tc - ac_tg - ac_tt) + tt_ac
    ag2 = tt_ac - ag_ac
    at2 = tt_ac - at_ac
    ca2 = tt_ac - ca_ac
    cc2 = tt_ac - cc_ac
    cg2 = tt_ac - cg_ac
    ct2 = tt_ac - ct_ac
    ga2 = tt_ac - ga_ac
    gc2 = tt_ac - gc_ac
    gg2 = tt_ac - gg_ac
    gt2 = tt_ac - gt_ac
    ta2 = tt_ac - ta_ac
    tc2 = tt_ac - tc_ac
    tg2 = tt_ac - tg_ac
    constant2 = tt_ac
    equation2 = [aa2, ac2, ag2, at2, ca2, cc2, cg2, ct2, ga2, gc2, gg2, gt2, ta2, tc2, tg2]

    #Equation 3
    aa3 = tt_ag - aa_ag
    ac3 = tt_ag - ac_ag
    ag3 = 1 - (1 - ag_aa - ag_ac - ag_at - ag_ca - ag_cc - ag_cg - ag_ct - ag_ga - ag_gc - ag_gg - ag_gt - ag_ta - ag_tc - ag_tg - ag_tt) + tt_ag
    at3 = tt_ag - at_ag
    ca3 = tt_ag - ca_ag
    cc3 = tt_ag - cc_ag
    cg3 = tt_ag - cg_ag
    ct3 = tt_ag - ct_ag
    ga3 = tt_ag - ga_ag
    gc3 = tt_ag - gc_ag
    gg3 = tt_ag - gg_ag
    gt3 = tt_ag - gt_ag
    ta3 = tt_ag - ta_ag
    tc3 = tt_ag - tc_ag
    tg3 = tt_ag - tg_ag
    constant3 = tt_ag
    equation3 = [aa3, ac3, ag3, at3, ca3, cc3, cg3, ct3, ga3, gc3, gg3, gt3, ta3, tc3, tg3]

    #Equation 4
    aa4 = tt_at - aa_at
    ac4 = tt_at - ac_at
    ag4 = tt_at - ag_at
    at4 = 1 - (1 - at_aa - at_ac - at_ag - at_ca - at_cc - at_cg - at_ct - at_ga - at_gc - at_gg - at_gt - at_ta - at_tc - at_tg - at_tt) + tt_at
    ca4 = tt_at - ca_at
    cc4 = tt_at - cc_at
    cg4 = tt_at - cg_at
    ct4 = tt_at - ct_at
    ga4 = tt_at - ga_at
    gc4 = tt_at - gc_at
    gg4 = tt_at - gg_at
    gt4 = tt_at - gt_at
    ta4 = tt_at - ta_at
    tc4 = tt_at - tc_at
    tg4 = tt_at - tg_at
    constant4 = tt_at
    equation4 = [aa4, ac4, ag4, at4, ca4, cc4, cg4, ct4, ga4, gc4, gg4, gt4, ta4, tc4, tg4]

    #Equation 5
    aa5 = tt_ca - aa_ca
    ac5 = tt_ca - ac_ca
    ag5 = tt_ca - ag_ca
    at5 = tt_ca - at_ca
    ca5 = 1 - (1 - ca_aa - ca_ac - ca_ag - ca_at - ca_cc - ca_cg - ca_ct - ca_ga - ca_gc - ca_gg - ca_gt - ca_ta - ca_tc - ca_tg - ca_tt) + tt_ca
    cc5 = tt_ca - cc_ca
    cg5 = tt_ca - cg_ca
    ct5 = tt_ca - ct_ca
    ga5 = tt_ca - ga_ca
    gc5 = tt_ca - gc_ca
    gg5 = tt_ca - gg_ca
    gt5 = tt_ca - gt_ca
    ta5 = tt_ca - ta_ca
    tc5 = tt_ca - tc_ca
    tg5 = tt_ca - tg_ca
    constant5 = tt_ca
    equation5 = [aa5, ac5, ag5, at5, ca5, cc5, cg5, ct5, ga5, gc5, gg5, gt5, ta5, tc5, tg5]

    #Equation 6
    aa6 = tt_cc - aa_cc
    ac6 = tt_cc - ac_cc
    ag6 = tt_cc - ag_cc
    at6 = tt_cc - at_cc
    ca6 = tt_cc - ca_cc
    cc6 = 1 - (1 - cc_aa - cc_ac - cc_ag - cc_at - cc_ca - cc_cg - cc_ct - cc_ga - cc_gc - cc_gg - cc_gt - cc_ta - cc_tc - cc_tg - cc_tt) + tt_cc
    cg6 = tt_cc - cg_cc
    ct6 = tt_cc - ct_cc
    ga6 = tt_cc - ga_cc
    gc6 = tt_cc - gc_cc
    gg6 = tt_cc - gg_cc
    gt6 = tt_cc - gt_cc
    ta6 = tt_cc - ta_cc
    tc6 = tt_cc - tc_cc
    tg6 = tt_cc - tg_cc
    constant6 = tt_cc
    equation6 = [aa6, ac6, ag6, at6, ca6, cc6, cg6, ct6, ga6, gc6, gg6, gt6, ta6, tc6, tg6]

    #Equation 7
    aa7 = tt_cg - aa_cg
    ac7 = tt_cg - ac_cg
    ag7 = tt_cg - ag_cg
    at7 = tt_cg - at_cg
    ca7 = tt_cg - ca_cg
    cc7 = tt_cg - cc_cg
    cg7 = 1 - (1 - cg_aa - cg_ac - cg_ag - cg_at - cg_ca - cg_cc - cg_ct - cg_ga - cg_gc - cg_gg - cg_gt - cg_ta - cg_tc - cg_tg - cg_tt) + tt_cg
    ct7 = tt_cg - ct_cg
    ga7 = tt_cg - ga_cg
    gc7 = tt_cg - gc_cg
    gg7 = tt_cg - gg_cg
    gt7 = tt_cg - gt_cg
    ta7 = tt_cg - ta_cg
    tc7 = tt_cg - tc_cg
    tg7 = tt_cg - tg_cg
    constant7 = tt_cg
    equation7 = [aa7, ac7, ag7, at7, ca7, cc7, cg7, ct7, ga7, gc7, gg7, gt7, ta7, tc7, tg7]

    #Equation 8
    aa8 = tt_ct - aa_ct
    ac8 = tt_ct - ac_ct
    ag8 = tt_ct - ag_ct
    at8 = tt_ct - at_ct
    ca8 = tt_ct - ca_ct
    cc8 = tt_ct - cc_ct
    cg8 = tt_ct - cg_ct
    ct8 = 1 - (1 - ct_aa - ct_ac - ct_ag - ct_at - ct_ca - ct_cc - ct_cg - ct_ga - ct_gc - ct_gg - ct_gt - ct_ta - ct_tc - ct_tg - ct_tt) + tt_ct
    ga8 = tt_ct - ga_ct
    gc8 = tt_ct - gc_ct
    gg8 = tt_ct - gg_ct
    gt8 = tt_ct - gt_ct
    ta8 = tt_ct - ta_ct
    tc8 = tt_ct - tc_ct
    tg8 = tt_ct - tg_ct
    constant8 = tt_ct
    equation8 = [aa8, ac8, ag8, at8, ca8, cc8, cg8, ct8, ga8, gc8, gg8, gt8, ta8, tc8, tg8]

    #Equation 9
    aa9 = tt_ga - aa_ga
    ac9 = tt_ga - ac_ga
    ag9 = tt_ga - ag_ga
    at9 = tt_ga - at_ga
    ca9 = tt_ga - ca_ga
    cc9 = tt_ga - cc_ga
    cg9 = tt_ga - cg_ga
    ct9 = tt_ga - ct_ga
    ga9 = 1 - (1 - ga_aa - ga_ac - ga_ag - ga_at - ga_ca - ga_cc - ga_cg - ga_ct - ga_gc - ga_gg - ga_gt - ga_ta - ga_tc - ga_tg - ga_tt) + tt_ga
    gc9 = tt_ga - gc_ga
    gg9 = tt_ga - gg_ga
    gt9 = tt_ga - gt_ga
    ta9 = tt_ga - ta_ga
    tc9 = tt_ga - tc_ga
    tg9 = tt_ga - tg_ga
    constant9 = tt_ga
    equation9 = [aa9, ac9, ag9, at9, ca9, cc9, cg9, ct9, ga9, gc9, gg9, gt9, ta9, tc9, tg9]

    #Equation 10
    aa10 = tt_gc - aa_gc
    ac10 = tt_gc - ac_gc
    ag10 = tt_gc - ag_gc
    at10 = tt_gc - at_gc
    ca10 = tt_gc - ca_gc
    cc10 = tt_gc - cc_gc
    cg10 = tt_gc - cg_gc
    ct10 = tt_gc - ct_gc
    ga10 = tt_gc - ga_gc
    gc10 = 1 - (1 - gc_aa - gc_ac - gc_ag - gc_at - gc_ca - gc_cc - gc_cg - gc_ct - gc_ga - gc_gg - gc_gt - gc_ta - gc_tc - gc_tg - gc_tt) + tt_gc
    gg10 = tt_gc - gg_gc
    gt10 = tt_gc - gt_gc
    ta10 = tt_gc - ta_gc
    tc10 = tt_gc - tc_gc
    tg10 = tt_gc - tg_gc
    constant10 = tt_gc
    equation10 = [aa10, ac10, ag10, at10, ca10, cc10, cg10, ct10, ga10, gc10, gg10, gt10, ta10, tc10, tg10]

    #Equation 11
    aa11 = tt_gg - aa_gg
    ac11 = tt_gg - ac_gg
    ag11 = tt_gg - ag_gg
    at11 = tt_gg - at_gg
    ca11 = tt_gg - ca_gg
    cc11 = tt_gg - cc_gg
    cg11 = tt_gg - cg_gg
    ct11 = tt_gg - ct_gg
    ga11 = tt_gg - ga_gg
    gc11 = tt_gg - gc_gg
    gg11 = 1 - (1 - gg_aa - gg_ac - gg_ag - gg_at - gg_ca - gg_cc - gg_cg - gg_ct - gg_ga - gg_gc - gg_gt - gg_ta - gg_tc - gg_tg - gg_tt) + tt_gg
    gt11 = tt_gg - gt_gg
    ta11 = tt_gg - ta_gg
    tc11 = tt_gg - tc_gg
    tg11 = tt_gg - tg_gg
    constant11 = tt_gg
    equation11 = [aa11, ac11, ag11, at11, ca11, cc11, cg11, ct11, ga11, gc11, gg11, gt11, ta11, tc11, tg11]

    #Equation 12
    aa12 = tt_gt - aa_gt
    ac12 = tt_gt - ac_gt
    ag12 = tt_gt - ag_gt
    at12 = tt_gt - at_gt
    ca12 = tt_gt - ca_gt
    cc12 = tt_gt - cc_gt
    cg12 = tt_gt - cg_gt
    ct12 = tt_gt - ct_gt
    ga12 = tt_gt - ga_gt
    gc12 = tt_gt - gc_gt
    gg12 = tt_gt - gg_gt
    gt12 = 1 - (1 - gt_aa - gt_ac - gt_ag - gt_at - gt_ca - gt_cc - gt_cg - gt_ct - gt_ga - gt_gc - gt_gg - gt_ta - gt_tc - gt_tg - gt_tt) + tt_gt
    ta12 = tt_gt - ta_gt
    tc12 = tt_gt - tc_gt
    tg12 = tt_gt - tg_gt
    constant12 = tt_gt
    equation12 = [aa12, ac12, ag12, at12, ca12, cc12, cg12, ct12, ga12, gc12, gg12, gt12, ta12, tc12, tg12]

    #Equation 13
    aa13 = tt_ta - aa_ta
    ac13 = tt_ta - ac_ta
    ag13 = tt_ta - ag_ta
    at13 = tt_ta - at_ta
    ca13 = tt_ta - ca_ta
    cc13 = tt_ta - cc_ta
    cg13 = tt_ta - cg_ta
    ct13 = tt_ta - ct_ta
    ga13 = tt_ta - ga_ta
    gc13 = tt_ta - gc_ta
    gg13 = tt_ta - gg_ta
    gt13 = tt_ta - gt_ta
    ta13 = 1 - (1 - ta_aa - ta_ac - ta_ag - ta_at - ta_ca - ta_cc - ta_cg - ta_ct - ta_ga - ta_gc - ta_gg - ta_gt - ta_tc - ta_tg - ta_tt) + tt_ta
    tc13 = tt_ta - tc_ta
    tg13 = tt_ta - tg_ta
    constant13 = tt_ta
    equation13 = [aa13, ac13, ag13, at13, ca13, cc13, cg13, ct13, ga13, gc13, gg13, gt13, ta13, tc13, tg13]

    #Equation 14
    aa14 = tt_tc - aa_tc
    ac14 = tt_tc - ac_tc
    ag14 = tt_tc - ag_tc
    at14 = tt_tc - at_tc
    ca14 = tt_tc - ca_tc
    cc14 = tt_tc - cc_tc
    cg14 = tt_tc - cg_tc
    ct14 = tt_tc - ct_tc
    ga14 = tt_tc - ga_tc
    gc14 = tt_tc - gc_tc
    gg14 = tt_tc - gg_tc
    gt14 = tt_tc - gt_tc
    ta14 = tt_tc - ta_tc
    tc14 = 1 - (1 - tc_aa - tc_ac - tc_ag - tc_at - tc_ca - tc_cc - tc_cg - tc_ct - tc_ga - tc_gc - tc_gg - tc_gt - tc_ta - tc_tg - tc_tt) + tt_tc
    tg14 = tt_tc - tg_tc
    constant14 = tt_tc
    equation14 = [aa14, ac14, ag14, at14, ca14, cc14, cg14, ct14, ga14, gc14, gg14, gt14, ta14, tc14, tg14]

    #Equation 15
    aa15 = tt_tg - aa_tg
    ac15 = tt_tg - ac_tg
    ag15 = tt_tg - ag_tg
    at15 = tt_tg - at_tg
    ca15 = tt_tg - ca_tg
    cc15 = tt_tg - cc_tg
    cg15 = tt_tg - cg_tg
    ct15 = tt_tg - ct_tg
    ga15 = tt_tg - ga_tg
    gc15 = tt_tg - gc_tg
    gg15 = tt_tg - gg_tg
    gt15 = tt_tg - gt_tg
    ta15 = tt_tg - ta_tg
    tc15 = tt_tg - tc_tg
    tg15 = 1 - (1 - tg_aa - tg_ac - tg_ag - tg_at - tg_ca - tg_cc - tg_cg - tg_ct - tg_ga - tg_gc - tg_gg - tg_gt - tg_ta - tg_tc - tg_tt) + tt_tc
    constant15 = tt_tg
    equation15 = [aa15, ac15, ag15, at15, ca15, cc15, cg15, ct15, ga15, gc15, gg15, gt15, ta15, tc15, tg15]

    #Create numpy arrays
    equations = np.array([equation1, equation2, equation3, equation4, equation5,
                    equation6, equation7, equation8, equation9, equation10, equation11,
                    equation12, equation13, equation14, equation15
                    ])

    answer = np.array([constant1, constant2, constant3, constant4, constant5,
                    constant6, constant7, constant8, constant9, constant10, constant11,
                    constant12, constant13, constant14, constant15])

    #Solve the simultaneous equations - note, we estimate tt equilibrium as 1 - the others
    solved = np.linalg.solve(equations,answer)
    aa_eq = float(solved[0])
    ac_eq = float(solved[1])
    ag_eq = float(solved[2])
    at_eq = float(solved[3])
    ca_eq = float(solved[4])
    cc_eq = float(solved[5])
    cg_eq = float(solved[6])
    ct_eq = float(solved[7])
    ga_eq = float(solved[8])
    gc_eq = float(solved[9])
    gg_eq = float(solved[10])
    gt_eq = float(solved[11])
    ta_eq = float(solved[12])
    tc_eq = float(solved[13])
    tg_eq = float(solved[14])
    tt_eq = 1 - aa_eq - ac_eq - ag_eq - at_eq - ca_eq - cc_eq - cg_eq - ct_eq - ga_eq - gc_eq - gg_eq - gt_eq - ta_eq - tc_eq - tg_eq

    #Estimate equilibrium GC content
    predicted_gc = (0.5 * (ac_eq + ag_eq + ca_eq + ct_eq + ga_eq + gt_eq + tc_eq + tg_eq)) + cc_eq + cg_eq + gc_eq + gg_eq

    return aa_eq, ac_eq, ag_eq, at_eq, ca_eq, cc_eq, cg_eq, ct_eq, ga_eq, gc_eq, gg_eq, gt_eq, ta_eq, tc_eq, tg_eq, tt_eq, predicted_gc

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h

### RUN ###

if __name__ == '__main__':
  main()
