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

# deviation_effort_fishing_mortality
D_i <- c(-0.23, 0.11, 0.11, -0.44, -0.31, 0.09, 0.11, -0.12, -0.35, -0.08, -0.09, -0.07, 0.20, 0.47, 0.30, 0.19, -0.20, 0.29, -0.46, -0.08)

# errors in total catch estimates <= 乱数発生
epsilon_i <- round(rnorm(mean=0, sd=0.1, 20), digits=2)
# 再解析しても結果が一致するように保存	
epsilon_i <- c(-0.14, -0.13, -0.03, -0.07, 0.07, -0.17, 0.05, 0.11, 0.01, -0.09, -0.05, -0.04, 0.13, -0.04, -0.02, 0.11, 0.10, 0.04, -0.08, -0.05)

# natural mortality rate
M_ij <- 0.5
# standard deviation for errors in total catch estimates
sigma <- 0.1
# standard deviation for deviation from average fishing mortality
sigma_1 <- 0.1
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



# calculaltion for beta_ij
# i 1:20, j 4:13 (data.frameでは1:10)
calc_beta_ij <- function(i, j, b_1, b_2, s, r) {
	O_i <- sum(CAA[i,])
	O_i1 <- sum(CAA[i+1,])
	j_s <- -1 + 2*(1 - s^(j-1))/(1 - s^(r-1))
	j_s_1 <- -1 + 2*(1 - s^(j))/(1 - s^(r-1))
	F_ij <- exp(b_1*j_s + b_2*(j_s^2) + D_i[i])
	F_i1j1 <- exp(b_1*j_s_1 + b_2*(j_s_1^2) + D_i[i+1])
	M_i1j1 <- M_ij
	beta_ij <- log(F_ij) - log(F_i1j1) - log(F_ij + M_ij) + log(F_i1j1 + M_i1j1) + 
	log(exp(F_ij + M_ij) - 1) - log(1 - exp(-F_i1j1 - M_i1j1)) - log(O_i) + log(O_i1) + 
	log(observed_P_ij[i+1,j-2]/exp(epsilon_i[i+1]))
	return(beta_ij)
}

# log-likelihood function
log_likelihood_1.9 <- function(x, y, z) {
	return(function(par) {
		b_1 <- par[1]
		b_2 <- par[2]
		s <- par[3]

		# sigma_exp_beta_ij data.frame作成
		for (i in 1:n-1) {
			sigma_exp_beta_j <- 0
			for (j in a:r-1) {
				beta_ij = calc_beta_ij(i=i, j=j, b_1=b_1, b_2=b_2, s=s, r=r)
				sigma_exp_beta_j <- sigma_exp_beta_j + exp(beta_ij)
			}
			if (i == 1) {
				sigma_exp_beta_ij <- sigma_exp_beta_j
			} else {
				sigma_exp_beta_ij <- append(sigma_exp_beta_ij, sigma_exp_beta_j)
			}
		}

		# likelihood
		LL1 <- 0
		for (i in 1:n-1) {
			for (j in a:r-1) {
				beta_ij <- calc_beta_ij(i=i, j=j, b_1=b_1, b_2=b_2, s=s, r=13)
				LL1 <- S_ij_400[i,j-3] * (beta_ij - log(sigma_exp_beta_ij[i]))
			}
		}

		LL2 <- 0
		for (i in 1:n-1) {
			LL2 <- 1/2 * log((sigma_exp_beta_ij[i]^2) / sigma^2)
		}

		LL3 <- m * log(sigma)

		LL <- LL1 - LL2 - LL3

		return(LL)
	})
}

# optimization
param_initial <- c(1.0, 1.0, 0.5)
res <- optim(par=param_initial, fn=log_likelihood_1.9(param_initial), method="SANN", hessian=TRUE, control=list(fnscale=-1))

# 推定したパラメーターを初期値にしてもう一度解析
param_initial <- c(-7.749316, 6.716972, 0.605138)
res2 <- optim(par=param_initial, fn=log_likelihood_1.9(param_initial), method="SANN", hessian=TRUE, control=list(fnscale=-1))

# パラメーターをセット
b_1 <- -7.749316
b_2 <- 6.716972
s <- 0.605138

# b_1, b_2, sからF_ijのdata.frame作成
for (i in 1:n) {
	for (j in a:r) {
		j_s <- -1 + 2*(1 - s^(j-1))/(1 - s^(r-1))
		F_j <- exp(b_1*j_s + b_2*(j_s^2) + D_i[i])

		# data.frameに入れる
		if (j == a) {
			F_i <- F_j
		} else {
			F_i <- append(F_i, F_j)
		}
	}
	if (i == 1) {
		F_ij <- t(data.frame(F_i))
		colnames(F_ij) <- c(a:r)
	} else {
		F_ij <- rbind(F_ij, F_i)
	}
}
rownames(F_ij) <- c(1:n)


# F_ij => NAA
M_ij <- 0.5
for (i in 1:n) {
	for (j in a:r) {
		Fij <- F_ij[i,(j-3)]
		C_ij <- CAA[i,(j-3)]
		# equation 1.3
		N_ij = (C_ij*(Fij + M_ij))/(Fij*(1 - exp(-Fij - M_ij)))

		if (j == a) {
			N_i <- N_ij
		} else {
			N_i <- append(N_i, N_ij)
		}
	}
	# Table of the catch at age
	if (i == 1) {
		NAA <- data.frame(round(N_i, digits=0))
		NAA <- t(NAA)
		colnames(NAA) <- c(a:r)
	} else {
		NAA <- rbind(NAA, round(N_i, digits=0))
	}
}
rownames(NAA) <- c(1:n)



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