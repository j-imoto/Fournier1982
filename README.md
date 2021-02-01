# Fournier&Archibald,1982の再現を目指す
 
まだ作成途中です。
論文を読んでも理解できていない部分があります。 
 
# SET1
 
Catch at age (Table 2)、Aging error probability (Table 3)、Observed number at age in the aging samples (Table 4)から漁獲死亡率、リッカーのαとδ、Number at ageを推定する
自然死亡率: M = 0.5
 
* T2_CAA.tsv
* T3_aging_error_probability.tsv
* T4_observed_number_at_age_400.tsv
* set1+random_error_M0.5.R
 
# SET1
 
Catch at age (Table 2)、Observed number at age in the aging samples (Table 4)から漁獲死亡率、リッカーのαとδ、Number at ageを推定する
自然死亡率: M = 0.5
 
* T2_CAA.tsv
* T4_observed_number_at_age_400.tsv
* set1+random_error_M0.5.R

# Calc_NAA2-20
 
1年目のNAA(Table 1: 1行目)と漁獲死亡率(Table 5)からNumber at ageを計算する
 
* T5_instantaneous_fishing_mortality.tsv
* Calc_NAA2-20.R
 
# General script
 
frasyr、SAMと比較できるようにgeneralなscriptに書き換え中です
