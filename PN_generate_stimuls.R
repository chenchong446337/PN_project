## generating random sequences of stimulus
n_test <- rep(c("Pain", "Ctrl"), each= 10)
ID_mouse <- c("m3", "m7", "m16", "m17", "m18")

dat_n_test <- NULL
for (i in seq_along(ID_mouse)) {
  x <- sample(n_test)
  dat_n_test <- cbind(dat_n_test, x)
  
}
colnames(dat_n_test) <- ID_mouse

write.csv(dat_n_test, "random_pinprick.csv")
