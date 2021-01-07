setwd("OneDrive - 国立研究開発法人 水産研究・教育機構/ref/資源管理/Fournier&Archibald,1982_R/")
#setwd("~/Fournier&Archibald,1982_R/")

# NAA
# Year 1 ベクトル  <=  Table 1 The actual number of fish at age.
yr1_NAA <- c(400, 1000, 600, 200, 300, 50, 10, 20, 5, 2)

# mortality  <=  Table 5
instantaneous_fishing_mortality <- read.table("T5_instantaneous_fishing_mortality.tsv", sep="\t", header=TRUE)
rownames(instantaneous_fishing_mortality) <- 1:20
colnames(instantaneous_fishing_mortality) <- 4:13

deviation_effort_fishing_mortality <- c(-0.23, 0.11, 0.11, -0.44, -0.31, 0.09, 0.11, -0.12, -0.35, -0.08, -0.09, -0.07, 0.20, 0.47, 0.30, 0.19, -0.20, 0.29, -0.46, -0.08)

M <- 0.5	# natural_mortality_rate

# recruitment  <=  Table 6 (Ricker)
alpha = 9.0
delta = 0.006
# relative_reprcdustive_potential <- c(0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.0, 1.0, 1.0, 1.0)
recruitment_2_4 <- c(1000, 700, 300)  #  <=  Table 1
predicted_recruitment_5_20 <- c(50, 105, 153, 257, 493, 549, 404, 432, 551, 543, 548, 546, 536, 537, 542, 491)
error_in_stock_recruitment_5_20 <- c(0.177, -0.086, -0.268, -0.066, 0.713, -0.600, 0.167, 0.082, -0.378, 0.523, -0.429, 0.455, 0.319, 0.158, -0.152, -0.059)


# 全年度用のデータフレーム 1年目だけで2年目以降はまだ入っていない
NAA <- t(data.frame(yr1_NAA))
colnames(NAA) <- c(4:13)

for (i in 1:20) {
	for (j in 4:13) {
		# instantaneous_fishing_mortalityを[i,(j-3)]を見やすくする
		F_ij <- instantaneous_fishing_mortality[i,(j-3)]
		# deviationを入れる
		# F_ij <- instantaneous_fishing_mortality[i,(j-3)] * exp(deviation_effort_fishing_mortality[i])

		# Catch at age  eq1.3
		CAA_ij <- (F_ij/(F_ij + M)) * (1 -exp(-F_ij -M)) * NAA[i,(j-3)]
		if (j == 4) {
			CAA_i <- CAA_ij
		} else {
			CAA_i <- append(CAA_i, CAA_ij)
		}

		# Number at age i+1, j+1  eq1.4
		N_i1j1 <- exp(-(F_ij) -M) * NAA[i,(j-3)]
		if (j == 4) {
			N_i1 <- N_i1j1
		} else {
			N_i1 <- append(N_i1, N_i1j1)
		}
	}
	# recruitment
	if (i < 4) {	#  <=  Table 1
		recruitment_i1 <- recruitment_2_4[i]
	} else {		#  <=  Table 5  predicted recruitment & error in recruitment
		recruitment_i1 <- predicted_recruitment_5_20[i-3] * exp(error_in_stock_recruitment_5_20[i-3])
	}

	# 新しい年i+1の行として追加
	if (i < 20) {
		year_i1 <- c(recruitment_i1, N_i1)
		NAA <- rbind(NAA, round(year_i1[-11], digits=0))
	}

	# Table of the catch at age
	if (i == 1) {
		CAA <- data.frame(round(CAA_i, digits=0))
		CAA <- t(CAA)
		colnames(CAA) <- c(4:13)
	} else {
		CAA <- rbind(CAA, round(CAA_i, digits=0))
	}

}
rownames(NAA) <- c(1:20)
rownames(CAA) <- c(1:20)
