####plots for distribution - plate notation


plot_dist <- function(dist, labels=c(), scale = 1, color="skyblue", plot_dist_name=T) {
  old_par <- par(mar = c(0.3, 0, 0, 0), xaxt='n', yaxt='n',ann=FALSE, bty="n", xpd=NA)  
  x <- dist$x
  y <- do.call(dist$ddist, c(list(x=x), dist$ddist_params))
  # To always anchor the plot at zero and give some extra top space if neccecary.
  plot(c(x[1:2], x), c(0,  max(y) / (1- dist$top_space), y), type="l", col="transparent")
  
  # only draw where the distribution is not zero
  points_to_NA <- filter(c(0, y, 0), filter=c(1,1, 1)) == 0
  points_to_NA <- points_to_NA[-c(1, length(points_to_NA))]
  y[points_to_NA] <- NA
  if("bar" %in% dist$plot_type) {
    lines(x, y, type="h", col=color, lwd=6, lend=1)
    # Using legend to draw a white transparent box behind the text
    if(plot_dist_name) {
      legend(grconvertX(dist$name_pos[1], from="npc"), grconvertY(dist$name_pos[2], from="npc"),
             dist$name, cex=1.5 * scale, xjust=0.5, yjust=0.5, bty="o", box.lwd = 0, box.col="transparent",
             bg=rgb(1,1, 1,0.5),x.intersp=-1, y.intersp=0 , text.col="transparent")
    }
  }
  if("line" %in% dist$plot_type) {
    lines(x, y, type="l", col=color, lwd=3 * scale)
  }
  lines(grconvertX(c(0.037, (1 - 0.037)), from="npc"), grconvertY(c(-0.02,-0.02), from="npc"), lwd=2 * scale)
  if(plot_dist_name) {
    text(grconvertX(dist$name_pos[1], from="npc"), grconvertY(dist$name_pos[2], from="npc"), dist$name, cex=1.5 * scale)
  }
  
  if(is.character(names(labels))) {
    for(label_name in names(labels)) {
      xpos <- dist$labels[[label_name]][1]
      ypos <- dist$labels[[label_name]][2]
      label <- labels[label_name]
      text(grconvertX(xpos, from="npc"), grconvertY(ypos, from="npc"), label, cex=2 * scale)
    }
  } else {
    for(i in seq_along(labels)) {
      xpos <- dist$labels[[i]][1]
      ypos <- dist$labels[[i]][2]
      label <- labels[i]
      text(grconvertX(xpos, from="npc"), grconvertY(ypos, from="npc"), label, cex=2)
    }
  }
  par(old_par)
}

dists <- list(
  mvn_normal = list(
    name = expression('~ MVN'[k]),
    name_pos = c(0.5, 0.1),
    plot_type = "line",
    x = seq(-3.3, 3.3, 0.01),
    top_space = 0,
    ddist = dnorm,
    ddist_params = list(mean=0, sd=1),
    labels = list(mean = c(0.5, 0.3), right_sd = c(0.80, 0.5), left_sd = c(0.10, 0.5))
  ),
  gamma = list(
    name = "~ Gamma",
    name_pos = c(0.3, 0.1),
    plot_type = "line",
    x = seq(0, 2, 0.01),
    top_space = 0,
    ddist = dgamma,
    ddist_params = list(shape=1.3, rate=2.5),
    labels = list(params = c(0.60, 0.5))
  ),
  inv_wishart = list(
    name = "~ Inv-Wishart",
    name_pos = c(0.42, 0.1),
    plot_type = "line",
    x = seq(0, 1.5, 1.5),
    top_space = 0,
    ddist = function(x, shape, scale) {scale^shape / gamma(shape) * x^(-shape-1)*exp(-scale/x)},
    ddist_params = list(shape=3, scale=1),
    labels = list(mean = c(0.35, 0.3), params = c(0.65, 0.5))
  ),
 
  log_normal = list(
    name = expression(~ log-MVN[k]),
    name_pos = c(0.48, 0.1),
    plot_type = "line",
    x = seq(0, 1.7, 0.01),
    top_space = 0,
    ddist = dlnorm,
    ddist_params = list(meanlog=-0.3, sdlog=0.4),
    labels = list(mean = c(0.4, 0.3), right_sd = c(0.80, 0.5), left_sd = c(0.20, 0.5))
  ),
  log_mvn = list(
    name = expression(~ log-MVN[k]),
    name_pos = c(0.48, 0.1),
    plot_type = "line",
    x = seq(0, 1.7, 0.01),
    top_space = 0,
    ddist = dlnorm,
    ddist_params = list(meanlog=-0.3, sdlog=0.4),
    labels = list(mean = c(0.45, 0.3), right_sd = c(0.80, 0.5), left_sd = c(0.10, 0.5))
  )
)

plot_dist_svg <- function(dist, labels=c(), fname="", color="skyblue", plot_dist_name=T) {
  if(fname == "") {
    fname = paste(gsub("\\W", "", gsub("\\s", "_", dist$name)), ".svg", sep="")
  }
  svg(fname, width=2.25, height=1.688, bg="transparent")
  plot_dist(dist, labels, color=color, plot_dist_name=plot_dist_name)
  dev.off()
}

plot_dist_png <- function(dist, labels=c(), fname="", color="skyblue", plot_dist_name=T) {
  if(fname == "") {
    fname = paste(gsub("\\W", "", gsub("\\s", "_", dist$name)), ".png", sep="")
  }
  png(fname, width=165, height=123, bg="transparent", res=72, )
  plot_dist(dist, labels, color=color, plot_dist_name=plot_dist_name)
  dev.off()
}

# Function that renders text as an image. Useful for constructing images of equations. 
# See ?plotmath for examples and documentation

plot_text_svg <- function(expr, fname) {
  svg(fname, bg="transparent")
  plot.new()
  text(0.5, 0.5, expr)
  dev.off()
}

plot_text_png <- function(expr, fname, pointsize=32, width=640, height=480 ) {
  png(fname, bg="transparent", width=width, height=height, pointsize=pointsize)
  plot.new()
  text(0.5, 0.5, expr)
  dev.off()
}
