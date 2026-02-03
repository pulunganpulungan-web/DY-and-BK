# ============================
# 0. Load Required Libraries
# ============================
library(devtools)
library(zoo)
library(ConnectednessApproach)
library(vars)
library(corrplot)
library(readxl)
library(frequencyConnectedness)
library(urca)

# ============================
# 1. Data Preparation
# ============================
View(DATA_3)
names(DATA_3)
class(DATA_3)

# Rename first column to "Date"
names(DATA_3)[1] <- "Date"
DATA_3$Date <- as.Date(DATA_3$Date)

# Set column names
colnames(DATA_3) <- c("Date", "US", "UE", "Japan", "UK", "China", "Indonesia", "OP")

# Convert to zoo object
df1 <- DATA_3[, -1]
zoodataDATA_3 <- zoo(df1, order.by = DATA_3$Date)
View(zoodataDATA_3)


# ============================
# 2. Statistical descriptive
# ============================
DATA_3 <- read_excel("DATA 3.xlsx")
names(DATA_3)

#mean
colMeans(DATA_3[ , -1], na.rm = TRUE)
colMeans(DATA_3[ , c("US", "EU", "Japan", "UK", "China", "Indonesia", "op")], na.rm = TRUE)

#standard deviation
sd(DATA_3$US, na.rm = TRUE)
apply(DATA_3[ , c("US", "EU", "Japan", "UK", "China", "Indonesia", "op")], 2, sd, na.rm = TRUE)

#median
install.packages("matrixStats")  
library(matrixStats)
colMedians(as.matrix(DATA_3[ , c("US", "EU", "Japan", "UK", "China", "Indonesia")]), na.rm = TRUE)

#Min Value
apply(DATA_3[ , c("US", "EU", "Japan", "UK", "China", "Indonesia")], 2, min, na.rm = TRUE)


#Max Value
apply(DATA_3[ , c("US", "EU", "Japan", "UK", "China", "Indonesia")], 2, max, na.rm = TRUE)


#Correlation
install.packages("corrplot")
library(corrplot)
cor(DATA_3[ , c("US", "EU", "Japan", "UK", "China", "Indonesia")], use = "pairwise.complete.obs", method = "pearson")
#cor_matrix
write.csv(cor_matrix, file = "correlation_output.csv")


# ===================================
# 3. ADF Stationarity Test (Level)
# ===================================
cat("=== UJI ADF PADA DATA LEVEL ===\n")

for (col in colnames(zoodataDATA_3)) {
  cat("\nUji ADF untuk:", col, "\n")
  print(summary(ur.df(zoodataDATA_3[, col], type = "trend", selectlags = "AIC")))
}

# ============================
# 4. Differencing & ADF Re-Test
# ============================
zoodiff1 <- diff(zoodataDATA_3)

cat("\n\n=== UJI ADF PADA DATA DIFFERENCE ===\n")
for (col in colnames(zoodiff1)) {
  cat("\nUji ADF untuk:", col, "\n")
  print(summary(ur.df(zoodiff1[, col], type = "drift", selectlags = "AIC")))
}

View(zoodiff1)

# ============================
# 5. Optimal Lag Selection for VAR
# ============================
cat("\n\n=== PILIHAN LAG VAR ===\n")
lag_selection <- VARselect(zoodiff1, lag.max = 10, type = "const")
print(lag_selection$selection)

# ============================
# 6. Time-Domain Connectedness Analysis
# ============================
dyc_DATA_3 <- ConnectednessApproach(
  zoodiff1,
  nlag = 1,
  nfore = 10,
  window.size = 36,
  model = "VAR",
  connectedness = "Time"
)

# Output Tables & Plots
dyc_DATA_3$TABLE
ConnectednessApproach::PlotNetwork(dyc_DATA_3)
ConnectednessApproach::PlotFROM(dyc_DATA_3)
ConnectednessApproach::PlotTO(dyc_DATA_3, ylim = c(0, 100))
ConnectednessApproach::PlotTCI(dyc_DATA_3, ylim = c(0, 100))
ConnectednessApproach::PlotNET(dyc_DATA_3)


# Extract TCI from  rolling window (Time-Domain)
tci_rolling <- dyc_DATA_3$TCI

# Detais of TCI
head(tci_rolling)
tci_rolling
View(tci_rolling)
write.csv(tci_rolling, "TCI_Rolling.csv")
plot(tci_rolling, type = "l", col = "blue", lwd = 2,
     ylab = "Total Connectedness Index (TCI)", xlab = "Time",
     main = "Rolling Window Total Connectedness Index")

tci_rolling <- as.data.frame(dyc_DATA_3$TCI)
print(tci_rolling)


##Net pairwise for Indonesia
var_names <- dimnames(dyc_DATA_3$CT)[[1]]
n_var <- length(var_names)
n_win <- dim(dyc_DATA_3$CT)[3]
net_pairwise_total <- list()
net_pairwise_mean  <- list()
for (i in 1:(n_var - 1)) {
  for (j in (i + 1):n_var) {
    from_i_to_j <- dyc_DATA_3$CT[i, j, ]
    from_j_to_i <- dyc_DATA_3$CT[j, i, ]
    net_val <- from_i_to_j - from_j_to_i
    
    pair_name <- paste(var_names[i], "↔", var_names[j])
    
    # Total (jumlah semua window)
    net_pairwise_total[[pair_name]] <- sum(net_val)
    
    # Mean (rata-rata semua window)
    net_pairwise_mean[[pair_name]]  <- mean(net_val)
  }
}

net_pairwise_total_df <- data.frame(
  Pair = names(net_pairwise_total),
  Total_Net = unlist(net_pairwise_total)
)

net_pairwise_mean_df <- data.frame(
  Pair = names(net_pairwise_mean),
  Mean_Net = unlist(net_pairwise_mean)
)
net_pairwise_summary <- merge(net_pairwise_total_df, net_pairwise_mean_df, by = "Pair")

net_pairwise_summary_ID <- subset(net_pairwise_summary, grepl("Indonesia", Pair))

library(ConnectednessApproach)
library(reshape2)
library(ggplot2)
net_pairwise_list <- list()
for (i in 1:(n_var - 1)) {
  for (j in (i + 1):n_var) {
    # cek apakah pasangan mengandung Indonesia
    if (var_names[i] == "Indonesia" | var_names[j] == "Indonesia") {
      from_i_to_j <- dyc_DATA_3$CT[i, j, ]
      from_j_to_i <- dyc_DATA_3$CT[j, i, ]
      net_val <- from_i_to_j - from_j_to_i
      
      net_pairwise_list[[paste0(var_names[i], "_vs_", var_names[j])]] <-
        data.frame(
          Window = 1:n_win,
          NetPairwise = net_val,
          Pair = paste(var_names[i], "↔", var_names[j])
        )
    }
  }
}

net_pairwise_df <- do.call(rbind, net_pairwise_list)
ggplot(net_pairwise_df, aes(x = Window, y = NetPairwise)) +
  geom_line(color = "blue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~ Pair, scales = "free_y") +
  labs(
    title = "Net Pairwise Spillover untuk Pasangan dengan Indonesia",
    y = "Net Spillover (i→j minus j→i)",
    x = "Window"
  ) +
  theme_minimal()

net_pairwise_df$Year <- 1996 + (net_pairwise_df$Window - 1) / 12
ggplot(net_pairwise_df, aes(x = Year, y = NetPairwise)) +
  geom_line(color = "blue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~ Pair, scales = "free_y") +
  labs(
    title = "Net Pairwise Spillover untuk Pasangan dengan Indonesia",
    y = "Net Spillover (i→j minus j→i)",
    x = "Tahun"
  ) +
  theme_minimal()



# ============================
# 7. Frequency-Domain Connectedness Analysis
# ============================
partitionFIX <- c(pi+0.00001, pi/4, pi/12, 0)

dyc_freq <- ConnectednessApproach(
  zoodiff1,
  nlag = 1,
  nfore = 10,
  window.size = 36,
  model = "VAR",
  connectedness = "Frequency",
  Connectedness_config = list(
    FrequencyConnectedness = list(
      partition = partitionFIX,
      generalized = TRUE,
      scenario = "ABS"
    )
  )
)

# Output Summary & Tables
summary(dyc_freq)
dyc_freq$TABLE

# Frequency-Domain Plots
ConnectednessApproach::PlotNetwork(dyc_freq)
ConnectednessApproach::PlotFROM(dyc_freq)
ConnectednessApproach::PlotTO(dyc_freq, ylim = c(0, 100))
ConnectednessApproach::PlotNPDC(dyc_freq)
ConnectednessApproach::PlotNET(dyc_freq)
ConnectednessApproach::PlotTCI(dyc_freq, ylim = c(0, 100))


names(dyc_freq$TABLES)


# Frequency-Specific Metrics
dyc_freq$TCI              
dyc_freq$NET             
dyc_freq$TO              
dyc_freq$FROM             

# =================================
# 8. ROBUSTNESS CHECKS
# =================================
# TVP-VAR-EJC DI R (Antonakakis, N., Chatziantoniou, I., & Gabauer, D. (2020))
tvp_ejc <- ConnectednessApproach(
  zoodataDATA_3,     
  nlag = 2,
  nfore = 10,
  window.size = NULL,       
  model = "TVP-VAR",
  connectedness = "Time"
)

# Extract 
tvp_ejc$TABLE
tvp_values <- tvp_ejc$TCI
print(tvp_values)           


