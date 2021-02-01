# CAA: O_ij  <=  Table 2
CAA <- read.table("T2_CAA.tsv", sep="\t", header=TRUE)
rownames(CAA) <- 1:20
colnames(CAA) <- 4:13

# Aging error probability  <=  Table 3
aging_error_probability <- read.table("T3_aging_error_probability.tsv", sep="\t", header=TRUE)
rownames(aging_error_probability) <- 4:13

# Observed number at age: S_ij  <=  Table 4
S_ij_400 <- read.table("T4_observed_number_at_age_400.tsv", sep="\t", header=TRUE)
rownames(S_ij_400) <- 1:20
colnames(S_ij_400) <- 4:13

# age4:13 Relative reproductive potentials  <=  Table 6
f_j <- c(0.10, 0.20, 0.30, 0,50, 0.70, 0.90, 1.00, 1.00, 1.00, 1.00)

# square of deviation_effort_fishing_mortality
D_i <- c(-0.23, 0.11, 0.11, -0.44, -0.31, 0.09, 0.11, -0.12, -0.35, -0.08, -0.09, -0.07, 0.20, 0.47, 0.30, 0.19, -0.20, 0.29, -0.46, -0.08)
# square of deviation effort fishing mortality <= 乱数発生
D_i <- round(rnorm(mean=0, sd=0.09, 20), digits=2)

# square of errors in total catch estimates
epsilon_i <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
# square of errors in total catch estimates <= 乱数発生
epsilon_i <- round(rnorm(mean=0, sd=0.0025, 20), digits=4)

# square of errors in stock-recruitment relationship
epsilon_i_2 <- c(0.177, -0.086, -0.268, -0.066, 0.713, -0.600, 0.167, 0.082, -0.378, 0.523, -0.429, 0.455, 0.319, 0.158, -0.152, -0.059)
# square of errors in stock-recruitment relationship <= 乱数発生
epsilon_i_2 <- round(rnorm(mean=0, sd=0.09, 20), digits=2)

# standard deviation for errors in total catch estimates
sigma <- 0.0025
# standard deviation for deviation from average fishing mortality
sigma_1 <- 0.09
# standard deviation for deviation for the stock-recruitment relationship
sigma_2 <- 0.09
# penalty weights for 
w_0 <- 2
# penalty weights for fishing mortality
w_1 <- 2
# penalty weights for stock-recruitment relationship
w_2 <- 2
# natural mortality rate
M_ij <- 0.5
# Ricker
#alpha <- 0.9
#delta <- 0.006

# 年数
n <- 20
# 最終年齢
r <- 13
# 加入年齢
a <- 4
# 年齢を調べた魚の数 S_ij
m <- 400
# P_ijが0にならないように
small_number <- 1.0E-10

# observed P_ij   0にsmall_number: 1.0E-10を代入する
for (i in 1:n) {
	O_i <- sum(CAA[i,])
	for (j in a:r) {
		if (j == a) {
			if (CAA[i,j-3] == 0) {
				observed_P_ij_df <- small_number
			} else {
				observed_P_ij_df <- (CAA[i,j-3])/O_i
			}
		} else {
			if (CAA[i,(j-3)] == 0) {
				observed_P_ij_df <- append(observed_P_ij_df, small_number)
			} else {
				observed_P_ij_df <- append(observed_P_ij_df, (CAA[i,j-3])/O_i)
			}
		}
	}
	if (i == 1) {
		observed_P_ij <- t(data.frame(observed_P_ij_df))
		colnames(observed_P_ij) <- c(a:r)
	} else {
		observed_P_ij <- rbind(observed_P_ij, observed_P_ij_df)
	}
}
rownames(observed_P_ij) <- c(1:n)

# a_jk   a_jk[j-3] 年齢jのときの年齢推定正解率
a_jk <- aging_error_probability[,3]


# O_ij  <=  Table 2
# observed P_ij  <=  Table 2
# S_ij  <=  Table 4
# a_jk  <=  Table 3
# F_ij  <=  Table 5


# calculaltion for beta_ij
# i 1:20, j 4:13 (data.frameでは1:10)
calc_beta_ij <- function(i, j, b_1, b_2, s, r) {
	O_i <- sum(CAA[i,])
	O_i1 <- sum(CAA[i+1,])
	j_s <- -1 + 2*(1 - s^(j-1))/(1 - s^(r-1))
	j_s_1 <- -1 + 2*(1 - s^(j))/(1 - s^(r-1))
	# F_ij <- exp(b_1*j_s + b_2*(j_s^2) + D_i[i])
	F_ij <- exp(b_1*j_s + b_2*(j_s^2))
	# F_i1j1 <- exp(b_1*j_s_1 + b_2*(j_s_1^2) + D_i[i+1])
	F_i1j1 <- exp(b_1*j_s_1 + b_2*(j_s_1^2))
	M_i1j1 <- M_ij
	beta_ij <- log(F_ij) - log(F_i1j1) - log(F_ij + M_ij) + log(F_i1j1 + M_i1j1) + 
	log(exp(F_ij + M_ij) - 1) - log(1 - exp(-F_i1j1 - M_i1j1)) - log(O_i) + log(O_i1) + 
	log(observed_P_ij[i+1,j-2]/exp(epsilon_i[i+1]))
	return(beta_ij)
}

# log-likelihood function
log_likelihood_6.3_5.3 <- function(v, w, x, y, z) {
	return(function(par) {
		b_1 <- par[1]
		b_2 <- par[2]
		s <- par[3]
		alpha <- par[4]
		delta <- par[5]

		# sigma_exp_beta_ijを計算 data.frameに入れる
		for (i in 1:19) {
			sigma_exp_beta_j <- 0
			for (j in 4:12) {
				beta_ij <- calc_beta_ij(i=i, j=j, b_1=b_1, b_2=b_2, s=s, r=13)
				sigma_exp_beta_j <- sigma_exp_beta_j + exp(beta_ij)
			}
			if (i == 1) {
				sigma_exp_beta_ij <- sigma_exp_beta_j
			} else {
				sigma_exp_beta_ij <- append(sigma_exp_beta_ij, sigma_exp_beta_j)
			}
		}

		# likelihood 6.3
		LL1 <- 0
		for (i in 1:19) {
			for (j in 4:12) {
				for (k in 4:12) {
					beta_ij <- calc_beta_ij(i=i, j=j, b_1=b_1, b_2=b_2, s=s, r=13)
					sigma_k <- a_jk[k-3]*(exp(beta_ij)/sigma_exp_beta_ij[i])
				}
				LL1 <- S_ij_400[i,j-3] * log(sigma_k)
			}
		}

		LL2 <- 0
		for (i in 1:19) {
			LL2 <- -1/2 * (sigma_exp_beta_ij[i]^2 / sigma)
			# LL2 <- -1/2 * log(sigma_exp_beta_ij[i]^2)
		}

		# LL3 <- -19 * log(sigma)

		# effort-fishing mortality relationship
		LL4 <- 0
		for (i in 1:19) {
			LL4 <- -1/2 * D_i[i]^2
		}

		# likelihood 5.3 recruitment
		# n_ijを計算　data.frameに格納
		for (i in 1:20) {
			for (j in 4:13) {
			# f_ijを計算する
				j_s <- -1 + 2*(1 - s^(j-1))/(1 - s^12)
				f_ij <- exp(b_1*j_s + b_2*(j_s^2) + D_i[i])
				# n_ijを計算する
				n_j <- (CAA[i,j-3]*(f_ij + M_ij))/(f_ij*(1 - exp(-f_ij - M_ij)))
				if (j == 4) {
					n_i <-  n_j
				} else {
					n_i <- append(n_i, n_j)
				}
			}
			if (i == 1) {
				n_ij <- t(data.frame(n_i))
				colnames(n_ij) <- c(4:13)
			} else {
				n_ij <- rbind(n_ij, n_i)
			}
		}
		rownames(n_ij) <- c(1:20)

		# P_iを計算　data.frameに格納
		for (i in 1:20) {
			for (j in 4:13) {
				if (j == 4) {
					p_i <- f_j[j-3] * n_ij[i,j-3]
				} else {
					p_i <- p_i + (f_j[j-3] * n_ij[i,j-3])
				}
			}
			if (i == 1) {
				P_i <- p_i
			} else {
				P_i <- append(P_i, p_i)
			}
		}
		LL5 <- 0
		for (i in 1:16) {
			N_ia_1 <- alpha * P_i[i] * exp(-delta*P_i[i]) * exp(epsilon_i_2[i])
			N_ia_1 <- n_ij[i+4,1]
			# LL5 <- -(1/2 * (log(N_ia_1/P_i[i]) - log(alpha) + delta*P_i[i])^2 / sigma_2)
			LL5 <- -1/2 * (log(N_ia_1/P_i[i]) - log(alpha) + delta*P_i[i])^2
		}

		LL <- LL1 + LL2 + w_1*LL4 + w_2*LL5
		LL <- LL1 + LL2 + LL3 + w_1*LL4 + w_2*LL5
		# LL <- LL1 + w_0*LL2 + LL3 + w_1*LL4 + w_2*LL5

		return(LL)
	})
}

# optimization
#param_initial <- c(0.0, 0.0, 0.5, 9.0, 0.006)
param_initial <- c(0.0, 0.0, 0.5, 0.0, 0.0)
set1_res_M0.5 <- optim(par=param_initial, fn=log_likelihood_6.3_5.3(param_initial), method="SANN", hessian=TRUE, control=list(fnscale=-1))

#if (0) {
# 推定したパラメーターを初期値にしてもう一度解析
#param_1st_run <- c(set1_res_M0.5_1st$par[1], set1_res_M0.5_1st$par[2], set1_res_M0.5_1st$par[3], set1_res_M0.5_1st$par[4], set1_res_M0.5_1st$par[5])
#set1_res_M0.5 <- optim(par=param_1st_run, fn=log_likelihood_6.3_5.3(param_1st_run), method="SANN", hessian=TRUE, control=list(fnscale=-1))

# パラメーターをセット
b_1 <- set1_res_M0.5$par[1]
b_2 <- set1_res_M0.5$par[2]
s <- set1_res_M0.5$par[3]
#}

# パラメーターをセット
#b_1 <- set1_res_M0.5_1st$par[1]
#b_2 <- set1_res_M0.5_1st$par[2]
#s <- set1_res_M0.5_1st$par[3]

# b_1, b_2, sからF_ijのdata.frame作成
for (i in 1:20) {
	for (j in 4:13) {
		j_s <- -1 + 2*(1 - s^(j-1))/(1 - s^(r-1))
		# F_j <- exp(b_1*j_s + b_2*(j_s^2))
		F_j <- exp(b_1*j_s + b_2*(j_s^2) + D_i[i])

		# data.frameに入れる
		if (j == 4) {
			F_i <- F_j
		} else {
			F_i <- append(F_i, F_j)
		}
	}
	if (i == 1) {
		F_ij <- t(data.frame(F_i))
		colnames(F_ij) <- c(4:13)
	} else {
		F_ij <- rbind(F_ij, F_i)
	}
}
rownames(F_ij) <- c(1:20)
write.table(F_ij, "set1+random_error_Mij0.5_Fij.tsv", sep="\t")


# F_ij => NAA
M_ij <- 0.5
for (i in 1:20) {
	for (j in 4:13) {
		Fij <- F_ij[i,(j-3)]
		C_ij <- CAA[i,(j-3)]
		# equation 1.3
		N_ij <- (C_ij*(Fij + M_ij))/(Fij*(1 - exp(-Fij - M_ij)))

		if (j == 4) {
			N_i <- N_ij
		} else {
			N_i <- append(N_i, N_ij)
		}
	}
	# Table of the catch at age
	if (i == 1) {
		NAA <- data.frame(round(N_i, digits=0))
		NAA <- t(NAA)
		colnames(NAA) <- c(4:13)
	} else {
		NAA <- rbind(NAA, round(N_i, digits=0))
	}
}
rownames(NAA) <- c(1:20)
write.table(NAA, "set1+random_error_Mij0.5_NAA.tsv", sep="\t")

# print
print(NAA)
print(F_ij)
print(D_i)
print(epsilon_i)
print(epsilon_i_2)
print(set1_res_M0.5)

# ファイル出力
sink("set1+random_error_Mij0.5.txt")
print("NAA")
print(NAA)
print("F_ij")
print(F_ij)
print("D_i")
print(D_i)
print("epsilon_i")
print(epsilon_i)
print("epsilon_i_2")
print(epsilon_i_2)
print("set1_res_M0.5")
print(set1_res_M0.5)
sink()



#### 根拠となる数式 ####
if(0){

# Fournier&Archibald 1982 漁獲死亡率
equation 2.2
j_s = -1 + 2*(1 - s^(j-1))/(1 - s^(r-1))

equation 2.1
F_ij = exp(b_1*j_s + b_2*(j_s)^2 + D_i)


O_i = C_i * exp(ε_i)

β_ij = log(P_ij*C_i/O_i)
γ_i = log(O_i)

equation 1.8
β_ij + γ_i - β_i+1,j+1 - γ_i+1 = Φ_ij
β_ij = Φ_ij - γ_i + γ_i+1 + β_i+1,j+1

右辺のβ_ijとγ_iを書き換える
β_ij = Φ_ij - log(O_i) + log(O_i+1) + log(P_i+1,j+1*C_i+1/O_i+1)

equation 1.5
Φ_ij = log(F_ij) - log(F_i+1,j+1) - log(F_ij + M_ij) + log(F_i+1,j+1 + M_i+1,j+1) + log[exp(F_ij + M_ij) - 1] - log[1 - exp(-F_i+1,j+1 - M_i+1,j+1)]

β_ij = log(F_ij) - log(F_i+1,j+1) - log(F_ij + M_ij) + log(F_i+1,j+1 + M_i+1,j+1) + log[exp(F_ij + M_ij) - 1] - log[1 - exp(-F_i+1,j+1 - M_i+1,j+1)] - log(O_i) + log(O_i+1) + log(P_i+1,j+1*C_i+1/O_i+1)

O_i <= C_i, observed_P_ij <= P_ij
β_ij = log(F_ij) - log(F_i+1,j+1) - log(F_ij + M_ij) + log(F_i+1,j+1 + M_i+1,j+1) + log[exp(F_ij + M_ij) - 1] - log[1 - exp(-F_i+1,j+1 - M_i+1,j+1)] - log(O_i) + log(O_i+1) + log(observed_P_i+1,j+1)

equation 1.9
LL <- 0
for (i in 1:20) {
	sigma_exp_beta_ij <- 0
	for (j in 4:13) {
		j_s = -1 + 2(1 - s^(j-1))(1 - s^(r-1))
		F_ij = exp(b_1*j_s + b_2*(j_s)^2 + D_i[i])

		sigma_exp_beta_ij <- sigma_exp_beta_ij + exp(beta_ij)
	}
	for (j in 4:13) {
		LL <- LL + S_ij * (beta_ij - log(sigma_exp_beta_ij))
	}
	LL <- LL - 1/2 * log(sigma_exp_beta_ij)^2 / sigma^2
}
LL <- LL -m * log(sigma)

}
########