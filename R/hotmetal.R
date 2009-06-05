hotmetal <- function(n=64) {
  orig <- c("#010000", "#0C0000", "#170000", "#210000", "#2C0000",
            "#360000", "#410000", "#4C0000", "#560000", "#610000",
            "#6C0000", "#760000", "#810000", "#8B0000", "#960000",
            "#A10000", "#AB0000", "#B60000", "#C10000", "#CB0000",
            "#D60000", "#E00000", "#EB0000", "#F60000", "#FF0100",
            "#FF0C00", "#FF1700", "#FF2100", "#FF2C00", "#FF3600",
            "#FF4100", "#FF4C00", "#FF5600", "#FF6100", "#FF6C00",
            "#FF7600", "#FF8100", "#FF8B00", "#FF9600", "#FFA100",
            "#FFAB00", "#FFB600", "#FFC100", "#FFCB00", "#FFD600",
            "#FFE000", "#FFEB00", "#FFF600", "#FFFF02", "#FFFF12",
            "#FFFF22", "#FFFF32", "#FFFF42", "#FFFF52", "#FFFF62",
            "#FFFF72", "#FFFF81", "#FFFF91", "#FFFFA1", "#FFFFB1",
            "#FFFFC1", "#FFFFD1", "#FFFFE1", "#FFFFF1")
  if(n == 64) 
    return(orig)
  rgb.hot <- t(col2rgb(orig))
  temp <- matrix(NA, ncol=3, nrow=n)
  x <- seq(0,1,,64)
  xg <- seq(0,1,,n)
  for(k in 1:3) {
    ## hold <- splint(x, rgb.hot[,k], xg)
    hold <- interpSpline(x, rgb.hot[,k])
    hold <- predict(hold, xg)$y
    hold[hold < 0] <- 0
    hold[hold > 255] <- 255
    temp[,k] <- round(hold)
  }
  rgb(temp[,1], temp[,2], temp[,3], maxColorValue=255)
}
