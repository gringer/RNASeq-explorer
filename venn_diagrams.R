#!/usr/bin/Rscript


draw.venn <- function(sets, add = FALSE, percent = FALSE, ...){
  ## sets: array of vectors containing elements to intersect
  ## ...: additional arguments passed through to draw.X.venn function
  ## return value: the intersection sizes for each area in the plot
  
  
{
  draw.triple.venn <-
    function (area1, area2, area3, n12, n23, n13, n123, category = rep("", 3), rotation = 1, reverse = FALSE, euler.d = TRUE, scaled = TRUE,
              lwd = rep(2, 3), lty = rep("solid", 3), col = rep("black", 3), fill = NULL, alpha = rep(0.5, 3),
              label.col = rep("black", 7), cex = rep(1, 7), fontface = rep("arial", 7),
              fontfamily = rep("arial", 7), cat.pos = c(-40, 40, 180), cat.dist = c(0.05, 0.05, 0.025),
              cat.col = rep("black", 3), cat.cex = rep(1, 3), cat.fontface = rep("arial", 3),
              cat.fontfamily = rep("arial", 3), cat.just = list(c(0.5, 1), c(0.5, 1), c(0.5, 0)),
              cat.default.pos = "outer", cat.prompts = FALSE, rotation.degree = 0,
              rotation.centre = c(0.5, 0.5), ind = TRUE, sep.dist = 0.05, offset = 0, ...) {
      if (length(category) == 1) {
        cat <- rep(category, 3)
      }
      else if (length(category) != 3) {
        stop("Unexpected parameter length for 'category'")
      }
      if (length(lwd) == 1) {
        lwd <- rep(lwd, 3)
      }
      else if (length(lwd) != 3) {
        stop("Unexpected parameter length for 'lwd'")
      }
      if (length(lty) == 1) {
        lty <- rep(lty, 3)
      }
      else if (length(lty) != 3) {
        stop("Unexpected parameter length for 'lty'")
      }
      if (length(col) == 1) {
        col <- rep(col, 3)
      }
      else if (length(col) != 3) {
        stop("Unexpected parameter length for 'col'")
      }
      if (length(label.col) == 1) {
        label.col <- rep(label.col, 7)
      }
      else if (length(label.col) != 7) {
        stop("Unexpected parameter length for 'label.col'")
      }
      if (length(cex) == 1) {
        cex <- rep(cex, 7)
      }
      else if (length(cex) != 7) {
        stop("Unexpected parameter length for 'cex'")
      }
      if (length(fontface) == 1) {
        fontface <- rep(fontface, 7)
      }
      else if (length(fontface) != 7) {
        stop("Unexpected parameter length for 'fontface'")
      }
      if (length(fontfamily) == 1) {
        fontfamily <- rep(fontfamily, 7)
      }
      else if (length(fontfamily) != 7) {
        stop("Unexpected parameter length for 'fontfamily'")
      }
      if (length(fill) == 1) {
        fill <- rep(fill, 3)
      }
      else if (length(fill) != 3 & length(fill) != 0) {
        stop("Unexpected parameter length for 'fill'")
      }
      if (length(alpha) == 1) {
        alpha <- rep(alpha, 3)
      }
      else if (length(alpha) != 3 & length(alpha) != 0) {
        stop("Unexpected parameter length for 'alpha'")
      }
      if (length(cat.pos) == 1) {
        cat.pos <- rep(cat.pos, 3)
      }
      else if (length(cat.pos) != 3) {
        stop("Unexpected parameter length for 'cat.pos'")
      }
      if (length(cat.dist) == 1) {
        cat.dist <- rep(cat.dist, 3)
      }
      else if (length(cat.dist) != 3) {
        stop("Unexpected parameter length for 'cat.dist'")
      }
      if (length(cat.col) == 1) {
        cat.col <- rep(cat.col, 3)
      }
      else if (length(cat.col) != 3) {
        stop("Unexpected parameter length for 'cat.col'")
      }
      if (length(cat.cex) == 1) {
        cat.cex <- rep(cat.cex, 3)
      }
      else if (length(cat.cex) != 3) {
        stop("Unexpected parameter length for 'cat.cex'")
      }
      if (length(cat.fontface) == 1) {
        cat.fontface <- rep(cat.fontface, 3)
      }
      else if (length(cat.fontface) != 3) {
        stop("Unexpected parameter length for 'cat.fontface'")
      }
      if (length(cat.fontfamily) == 1) {
        cat.fontfamily <- rep(cat.fontfamily, 3)
      }
      else if (length(cat.fontfamily) != 3) {
        stop("Unexpected parameter length for 'cat.fontfamily'")
      }
      if (!(class(cat.just) == "list" & length(cat.just) == 3)) {
        stop("Unexpected parameter format for 'cat.just'")
      }
      else if (!(length(cat.just[[1]]) == 2 & length(cat.just[[2]]) ==
                   2 & length(cat.just[[3]]) == 2)) {
        stop("Unexpected parameter format for 'cat.just'")
      }
      if (euler.d == FALSE & scaled == TRUE) {
        stop("Uninterpretable parameter combination\nPlease set both euler.d = FALSE and scaled = FALSE to force Venn diagrams.")
      }
      if (offset > 1 | offset < 0) {
        stop("'Offset' must be between 0 and 1.  Try using 'rotation.degree = 180' to achieve offsets in the opposite direction.")
      }
      cat.pos <- cat.pos + rotation.degree
      a1 <- area1 - n12 - n13 + n123
      a2 <- n12 - n123
      a3 <- area2 - n12 - n23 + n123
      a4 <- n13 - n123
      a5 <- n123
      a6 <- n23 - n123
      a7 <- area3 - n13 - n23 + n123
      areas <- c(a1, a2, a3, a4, a5, a6, a7)
      if (euler.d) {
        special.code <- VennDiagram::decide.special.case(areas)
        function.name <- paste("draw.", special.code, sep = "")
        if (function.name %in% ls("package:VennDiagram")) {
          f1 <- get(function.name)
          rst <- f1(a1 = areas[1], a2 = areas[2], a3 = areas[3],
                    a4 = areas[4], a5 = areas[5], a6 = areas[6],
                    a7 = areas[7], category = category, reverse = reverse,
                    cat.default.pos = cat.default.pos, lwd = lwd,
                    lty = lty, col = col, label.col = label.col,
                    cex = cex, fontface = fontface, fontfamily = fontfamily,
                    cat.pos = cat.pos, cat.dist = cat.dist, cat.col = cat.col,
                    cat.cex = cat.cex, cat.fontface = cat.fontface,
                    cat.fontfamily = cat.fontfamily, cat.just = cat.just,
                    cat.prompts = cat.prompts, fill = fill, alpha = alpha,
                    ...)
          rst <- VennDiagram::adjust.venn(VennDiagram::rotate.venn.degrees(gList1 = rst,
                                                                           angle = rotation.degree, x.centre = rotation.centre[1],
                                                                           y.centre = rotation.centre[2]), ...)
          if (ind) {
            grid.draw(rst)
          }
          return(rst)
        }
      }
      rotated <- VennDiagram::rotate(areas, category, lwd, lty,
                                     col, label.col, cex, fontface, fontfamily, cat.col, cat.cex,
                                     cat.fontface, cat.fontfamily, alpha, rotation, reverse,
                                     fill)
      for (i in 1:length(areas)) {
        areas[i] <- rotated[[1]][i]
      }
      category <- rotated[[2]]
      lwd <- rotated$lwd
      lty <- rotated$lty
      col <- rotated$col
      label.col <- rotated$label.col
      cex <- rotated$cex
      fontface <- rotated$fontface
      fontfamily <- rotated$fontfamily
      cat.col <- rotated$cat.col
      cat.cex <- rotated$cat.cex
      cat.fontface <- rotated$cat.fontface
      cat.fontfamily <- rotated$cat.fontfamily
      fill <- rotated$fill
      alpha <- rotated$alpha
      areas.error <- c("a1 <- area1 - n12 - n13 + n123", "a2 <- n12 - n123",
                       "a3 <- area2 - n12 - n23 + n123", "a4 <- n13 - n123",
                       "a5 <- n123", "a6 <- n23 - n123", "a7 <- area3 - n13 - n23 + n123")
      for (i in 1:length(areas)) {
        if (areas[i] < 0) {
          stop(paste("Impossible:", areas.error[i], "produces negative area"))
        }
      }
      for (i in 1:length(areas)) {
        if (areas[i]) {
          scaled <- FALSE
        }
      }
      is.defaults <- TRUE
      if (is.expression(category)) {
        is.defaults <- FALSE
      }
      if (all(cat.default.pos != "outer", cat.default.pos != "text",
              !is.defaults, cat.prompts)) {
        print("No default location recognized.  Automatically changing to 'outer'")
        cat.default.pos <- "outer"
      }
      if (all(cat.default.pos == "outer", !is.defaults, cat.prompts)) {
        print("Placing category labels at default outer locations.  Use 'cat.pos' and 'cat.dist' to modify location.")
        print(paste("Current 'cat.pos':", cat.pos[1], "degrees,",
                    cat.pos[2], "degrees"))
        print(paste("Current 'cat.dist':", cat.dist[1], ",",
                    cat.dist[2]))
      }
      if (all(cat.default.pos == "text", !is.defaults, cat.prompts)) {
        print("Placing category labels at default text locations.  Use 'cat.pos' and 'cat.dist' to modify location.")
        print(paste("Current 'cat.pos':", cat.pos[1], "degrees,",
                    cat.pos[2], "degrees"))
        print(paste("Current 'cat.dist':", cat.dist[1], ",",
                    cat.dist[2]))
      }
      grob.list <- gList()
      if (!exists("overrideTriple")) {
        r1 <- sqrt(100/pi)
        r2 <- r1
        r3 <- r1
      }
      else { ##scales <- SAM
        r1 <- sqrt(area1/pi)
        r2 <- sqrt(area2/pi)
        r3 <- sqrt(area3/pi)
      }
      max.circle.size = max(c(r1,r2,r3))
      shrink.factor <- max.circle.size/r1
      r1 <- r1 * shrink.factor
      r2 <- r2 * shrink.factor
      r3 <- r3 * shrink.factor
      if (!exists("overrideTriple")) {
        a <- find.dist(100, 100, 40) * shrink.factor
        b <- a
        c <- a
      }
      else {
        a <- find.dist(area1, area2, n12) * shrink.factor
        b <- find.dist(area2, area3, n23) * shrink.factor
        c <- find.dist(area1, area3, n13) * shrink.factor
      }
      x.centres <- vector(mode = "numeric", length = 3)
      y.centres <- vector(mode = "numeric", length = 3)
      beta <- (a^2 + c^2 - b^2)/(2 * a * c)
      gamma <- sqrt(1 - beta^2)
      x.centres[1] <- (r1 - r2 - a + 1)/2
      x.centres[3] <- x.centres[1] + c * beta
      y.centres[3] <- (r3 - r1 + 1 - c * gamma)/2
      y.centres[1] <- y.centres[3] + c * gamma
      x.centres[2] <- x.centres[1] + a
      y.centres[2] <- y.centres[1]
      radii <- c(r1, r2, r3)
      for (i in 1:3) {
        grob.list <- gList(grob.list, VennDiagram::ellipse(x = x.centres[i], y = y.centres[i], a = radii[i], b = radii[i],
                                                           gp = gpar(lty = 0, fill = fill[i], alpha = alpha[i])))
      }
      for (i in 1:3) {
        grob.list <- gList(grob.list, VennDiagram::ellipse(x = x.centres[i], y = y.centres[i], a = radii[i], b = radii[i],
                                                           gp = gpar(lwd = lwd[i], lty = lty[i], col = col[i], fill = "transparent")))
      }
      new.x.centres <- vector(mode = "numeric", length = 3)
      new.y.centres <- vector(mode = "numeric", length = 3)
      cell.labels <- round(areas,2)
      cell.x <- vector(mode = "numeric", length = 7)
      cell.y <- vector(mode = "numeric", length = 7)
      x.cept.12 <- (r1^2 - r2^2 - x.centres[1]^2 + x.centres[2]^2)/(2 * (x.centres[2] - x.centres[1]))
      y.cept.12.1 <- sqrt(r1^2 - (x.cept.12 - x.centres[1])^2) + y.centres[1]
      y.cept.12.2 <- -sqrt(r1^2 - (x.cept.12 - x.centres[1])^2) + y.centres[1]
      theta <- acos((a^2 + c^2 - b^2)/(2 * a * c))
      new.x.centres[3] <- x.centres[1] + c
      l.x.cept.13 <- (r1^2 - r3^2 - x.centres[1]^2 + new.x.centres[3]^2)/(2 * (new.x.centres[3] - x.centres[1]))
      l.y.cept.13.1 <- sqrt(r1^2 - (l.x.cept.13 - x.centres[1])^2) + y.centres[1]
      l.y.cept.13.2 <- -sqrt(r1^2 - (l.x.cept.13 - x.centres[1])^2) + y.centres[1]
      rot <- sqrt(2 * r1^2 - 2 * r1^2 * cos(theta))
      x.cept.13.1 <- l.x.cept.13 + rot *
        cos(pi/2 - atan((l.y.cept.13.1 - y.centres[1])/(l.x.cept.13 - x.centres[1])) + theta/2)
      x.cept.13.2 <- l.x.cept.13 + rot *
        cos(pi/2 - atan((l.y.cept.13.2 - y.centres[1])/(l.x.cept.13 - x.centres[1])) + theta/2)
      y.cept.13.1 <- l.y.cept.13.1 - rot *
        sin(pi/2 - atan((l.y.cept.13.1 - y.centres[1])/(l.x.cept.13 - x.centres[1])) + theta/2)
      y.cept.13.2 <- l.y.cept.13.2 - rot *
        sin(pi/2 - atan((l.y.cept.13.2 - y.centres[1])/(l.x.cept.13 - x.centres[1])) + theta/2)
      theta <- -acos((a^2 + b^2 - c^2)/(2 * a * b))
      new.x.centres[3] <- x.centres[2] - b
      l.x.cept.23 <- (r2^2 - r3^2 - x.centres[2]^2 + new.x.centres[3]^2)/
        (2 * (new.x.centres[3] - x.centres[2]))
      l.y.cept.23.1 <- sqrt(r2^2 - (l.x.cept.23 - x.centres[2])^2) + y.centres[2]
      l.y.cept.23.2 <- -sqrt(r2^2 - (l.x.cept.23 - x.centres[2])^2) + y.centres[2]
      rot <- sqrt(2 * r2^2 - 2 * r2^2 * cos(theta))
      x.cept.23.1 <- l.x.cept.23 + rot * cos(pi/2 - atan((y.centres[2] - l.y.cept.23.1)/(x.centres[2] - l.x.cept.23)) + theta/2)
      x.cept.23.2 <- l.x.cept.23 + rot * cos(pi/2 - atan((y.centres[2] - l.y.cept.23.2)/(x.centres[2] - l.x.cept.23)) + theta/2)
      y.cept.23.1 <- l.y.cept.23.1 - rot * sin(pi/2 - atan((y.centres[2] - l.y.cept.23.1)/(x.centres[2] - l.x.cept.23)) + theta/2)
      y.cept.23.2 <- l.y.cept.23.2 - rot * sin(pi/2 - atan((y.centres[2] - l.y.cept.23.2)/(x.centres[2] - l.x.cept.23)) + theta/2)
      m <- (y.cept.23.2 - y.cept.23.1)/(x.cept.23.2 - x.cept.23.1)
      y.sect <- m * (x.cept.12 - x.cept.23.1) + y.cept.23.1
      cell.x[5] <- x.cept.12
      cell.y[5] <- y.sect
      m <- (y.cept.13.2 - y.cept.13.1)/(x.cept.13.2 - x.cept.13.1)
      y0 <- y.centres[2]
      x0 <- x.centres[2]
      b <- y.cept.13.1 - m * x.cept.13.1
      x.sect <- (m * y0 + x0 - m * b)/(m^2 + 1) + sqrt(r2^2 - ((y0 - m * x0 - b)/sqrt(1 + m^2))^2)/sqrt(1 + m^2)
      y.sect <- (m^2 * y0 + m * x0 + b)/(m^2 + 1) + m * sqrt(r2^2 - ((y0 - m * x0 - b)/sqrt(1 + m^2))^2)/sqrt(1 + m^2)
      cell.x[3] <- (x.cept.13.1 + x.sect)/2
      cell.y[3] <- (y.cept.13.1 + y.sect)/2
      m <- (y.cept.23.2 - y.cept.23.1)/(x.cept.23.2 - x.cept.23.1)
      y0 <- y.centres[1]
      x0 <- x.centres[1]
      b <- y.cept.23.1 - m * x.cept.23.1
      x.sect <- (m * y0 + x0 - m * b)/(m^2 + 1) - sqrt(r1^2 - ((y0 - m * x0 - b)/sqrt(1 + m^2))^2)/sqrt(1 + m^2)
      y.sect <- (m^2 * y0 + m * x0 + b)/(m^2 + 1) - m * sqrt(r1^2 - ((y0 - m * x0 - b)/sqrt(1 + m^2))^2)/sqrt(1 + m^2)
      cell.x[1] <- (x.cept.23.1 + x.sect)/2
      cell.y[1] <- (y.cept.23.1 + y.sect)/2
      y.sect <- -sqrt(r3^2 - (x.cept.12 - x.centres[3])^2) + y.centres[3]
      cell.x[7] <- x.cept.12
      cell.y[7] <- (y.cept.12.2 + y.sect)/2
      m <- (y.cept.23.2 - y.cept.23.1)/(x.cept.23.2 - x.cept.23.1)
      y0 <- y.centres[1]
      x0 <- x.centres[1]
      b <- y.cept.23.1 - m * x.cept.23.1
      x.sect <- (m * y0 + x0 - m * b)/(m^2 + 1) + sqrt(r1^2 - ((y0 - m * x0 - b)/sqrt(1 + m^2))^2)/sqrt(1 + m^2)
      y.sect <- (m^2 * y0 + m * x0 + b)/(m^2 + 1) + m * sqrt(r1^2 - ((y0 - m * x0 - b)/sqrt(1 + m^2))^2)/sqrt(1 + m^2)
      cell.x[6] <- (x.cept.23.2 + x.sect)/2
      cell.y[6] <- (y.cept.23.2 + y.sect)/2
      m <- (y.cept.13.2 - y.cept.13.1)/(x.cept.13.2 - x.cept.13.1)
      y0 <- y.centres[2]
      x0 <- x.centres[2]
      b <- y.cept.13.1 - m * x.cept.13.1
      x.sect <- (m * y0 + x0 - m * b)/(m^2 + 1) - sqrt(r2^2 - ((y0 - m * x0 - b)/sqrt(1 + m^2))^2)/sqrt(1 + m^2)
      y.sect <- (m^2 * y0 + m * x0 + b)/(m^2 + 1) - m * sqrt(r2^2 - ((y0 - m * x0 - b)/sqrt(1 + m^2))^2)/sqrt(1 + m^2)
      cell.x[4] <- (x.cept.13.2 + x.sect)/2
      cell.y[4] <- (y.cept.13.2 + y.sect)/2
      y.sect <- sqrt(r3^2 - (x.cept.12 - x.centres[3])^2) + y.centres[3]
      cell.x[2] <- x.cept.12
      cell.y[2] <- (y.cept.12.1 + y.sect)/2
      for (i in 1:7) {
        grob.list <-
          gList(grob.list, textGrob(label = cell.labels[i], x = cell.x[i], y = cell.y[i],
                                    gp = gpar(col = label.col[i], cex = cex[i], fontface = fontface[i], fontfamily = fontfamily[i])))
      }
      text.location.mapping <- c(1, 3, 7)
      for (i in 1:3) {
        if ("outer" == cat.default.pos) {
          this.cat.pos <- find.cat.pos(x = x.centres[i], y = y.centres[i],
                                       pos = cat.pos[i], dist = cat.dist[i], r = radii[i])
        }
        else if ("text" == cat.default.pos) {
          this.cat.pos <- find.cat.pos(x = cell.x[text.location.mapping[i]],
                                       y = cell.y[text.location.mapping[i]], pos = cat.pos[i],
                                       dist = cat.dist[i])
        }
        else {
          stop("Invalid setting of cat.default.pos")
        }
        grob.list <- gList(grob.list, textGrob(label = category[i],
                                               x = this.cat.pos$x, y = this.cat.pos$y, just = cat.just[[i]],
                                               gp = gpar(col = cat.col[i], cex = cat.cex[i], fontface = cat.fontface[i],
                                                         fontfamily = cat.fontfamily[i])))
      }
      grob.list <- VennDiagram::adjust.venn(VennDiagram::rotate.venn.degrees(gList1 = grob.list,
                                                                             angle = rotation.degree, x.centre = rotation.centre[1],
                                                                             y.centre = rotation.centre[2]), ...)
      if (ind) {
        grid.draw(grob.list)
      }
      return(grob.list)
    }
}

draw.quad.venn <- function (area1, area2, area3, area4, n12, n13, n14, n23, n24,
                            n34, n123, n124, n134, n234, n1234, category = rep("", 4),
                            lwd = rep(2, 4), lty = rep("solid", 4), col = rep("black", 4),
                            fill = NULL, alpha = rep(0.5, 4), label.col = rep("black", 15),
                            cex = rep(1, 15), fontface = rep("plain", 15), fontfamily = rep("serif", 15),
                            cat.pos = c(-15, 15, 0, 0), cat.dist = c(0.22, 0.22, 0.11, 0.11),
                            cat.col = rep("black", 4), cat.cex = rep(1, 4),
                            cat.fontface = rep("plain", 4), cat.fontfamily = rep("serif", 4),
                            cat.just = rep(list(c(0.5, 0.5)), 4), rotation.degree = 0,
                            rotation.centre = c(0.5, 0.5), ind = TRUE, ...)
{
  if (length(category) == 1) {
    cat <- rep(category, 4)
  }
  else if (length(category) != 4) {
    stop("Unexpected parameter length for 'category'")
  }
  if (length(lwd) == 1) {
    lwd <- rep(lwd, 4)
  }
  else if (length(lwd) != 4) {
    stop("Unexpected parameter length for 'lwd'")
  }
  if (length(lty) == 1) {
    lty <- rep(lty, 4)
  }
  else if (length(lty) != 4) {
    stop("Unexpected parameter length for 'lty'")
  }
  if (length(col) == 1) {
    col <- rep(col, 4)
  }
  else if (length(col) != 4) {
    stop("Unexpected parameter length for 'col'")
  }
  if (length(label.col) == 1) {
    label.col <- rep(label.col, 15)
  }
  else if (length(label.col) != 15) {
    stop("Unexpected parameter length for 'label.col'")
  }
  if (length(cex) == 1) {
    cex <- rep(cex, 15)
  }
  else if (length(cex) != 15) {
    stop("Unexpected parameter length for 'cex'")
  }
  if (length(fontface) == 1) {
    fontface <- rep(fontface, 15)
  }
  else if (length(fontface) != 15) {
    stop("Unexpected parameter length for 'fontface'")
  }
  if (length(fontfamily) == 1) {
    fontfamily <- rep(fontfamily, 15)
  }
  else if (length(fontfamily) != 15) {
    stop("Unexpected parameter length for 'fontfamily'")
  }
  if (length(fill) == 1) {
    fill <- rep(fill, 4)
  }
  else if (length(fill) != 4 & length(fill) != 0) {
    stop("Unexpected parameter length for 'fill'")
  }
  if (length(alpha) == 1) {
    alpha <- rep(alpha, 4)
  }
  else if (length(alpha) != 4 & length(alpha) != 0) {
    stop("Unexpected parameter length for 'alpha'")
  }
  if (length(cat.pos) == 1) {
    cat.pos <- rep(cat.pos, 4)
  }
  else if (length(cat.pos) != 4) {
    stop("Unexpected parameter length for 'cat.pos'")
  }
  if (length(cat.dist) == 1) {
    cat.dist <- rep(cat.dist, 4)
  }
  else if (length(cat.dist) != 4) {
    stop("Unexpected parameter length for 'cat.dist'")
  }
  if (length(cat.col) == 1) {
    cat.col <- rep(cat.col, 4)
  }
  else if (length(cat.col) != 4) {
    stop("Unexpected parameter length for 'cat.col'")
  }
  if (length(cat.cex) == 1) {
    cat.cex <- rep(cat.cex, 4)
  }
  else if (length(cat.cex) != 4) {
    stop("Unexpected parameter length for 'cat.cex'")
  }
  if (length(cat.fontface) == 1) {
    cat.fontface <- rep(cat.fontface, 4)
  }
  else if (length(cat.fontface) != 4) {
    stop("Unexpected parameter length for 'cat.fontface'")
  }
  if (length(cat.fontfamily) == 1) {
    cat.fontfamily <- rep(cat.fontfamily, 4)
  }
  else if (length(cat.fontfamily) != 4) {
    stop("Unexpected parameter length for 'cat.fontfamily'")
  }
  if (!(class(cat.just) == "list" & length(cat.just) == 4 & 
          length(cat.just[[1]]) == 2 & length(cat.just[[2]]) == 
          2 & length(cat.just[[3]]) == 2 & length(cat.just[[4]]) == 
          2)) {
    stop("Unexpected parameter format for 'cat.just'")
  }
  cat.pos <- cat.pos + rotation.degree
  a6 <- n1234
  a12 <- n123 - a6
  a11 <- n124 - a6
  a5 <- n134 - a6
  a7 <- n234 - a6
  a15 <- n12 - a6 - a11 - a12
  a4 <- n13 - a6 - a5 - a12
  a10 <- n14 - a6 - a5 - a11
  a13 <- n23 - a6 - a7 - a12
  a8 <- n24 - a6 - a7 - a11
  a2 <- n34 - a6 - a5 - a7
  a9 <- area1 - a4 - a5 - a6 - a10 - a11 - a12 - a15
  a14 <- area2 - a6 - a7 - a8 - a11 - a12 - a13 - a15
  a1 <- area3 - a2 - a4 - a5 - a6 - a7 - a12 - a13
  a3 <- area4 - a2 - a5 - a6 - a7 - a8 - a10 - a11
  areas <- c(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14)
  areas.error <- c("a1  <- area3 - a2 - a4 - a5 - a6 - a7 - a12 - a13",
                   "a2  <- n34 - a6 - a5 - a7", "a3  <- area4 - a2 - a5 - a6 - a7 - a8 - a10 - a11",
                   "a4  <- n13 - a6 - a5 - a12", "a5  <- n134 - a6", "a6  <- n1234",
                   "a7  <- n234 - a6", "a8  <- n24 - a6 - a7 - a11", "a9  <- area1 - a4 - a5 - a6 - a10 - a11 - a12 - a15",
                   "a10 <- n14 - a6 - a5 - a11", "a11 <- n124 - a6", "a12 <- n123 - a6",
                   "a15 <- n12 - a6 - a11 - a12", "a13 <- n23 - a6 - a7 - a12",
                   "a14 <- area2 - a6 - a7 - a8 - a11 - a12 - a13 - a15")
  for (i in 1:length(areas)) {
    if (areas[i] < 0) {
      stop(paste("Impossible:", areas.error[i], "produces negative area"))
    }
  }
  grob.list <- gList()
  ellipse.positions <- matrix(nrow = 4, ncol = 7)
  colnames(ellipse.positions) <- c("x", "y", "a", "b", "rotation",
                                   "fill.mapping", "line.mapping")
  ellipse.positions[1, ] <- c(0.65, 0.47, 0.35, 0.2, 45, 2, 4)
  ellipse.positions[2, ] <- c(0.35, 0.47, 0.35, 0.2, 135, 1, 1)
  ellipse.positions[3, ] <- c(0.5, 0.57, 0.33, 0.15, 45, 4, 3)
  ellipse.positions[4, ] <- c(0.5, 0.57, 0.35, 0.15, 135, 3, 2)
  for (i in 1:4) {
    grob.list <- gList(grob.list,
                       VennDiagram::ellipse(x = ellipse.positions[i, "x"],
                                            y = ellipse.positions[i, "y"],
                                            a = ellipse.positions[i, "a"],
                                            b = ellipse.positions[i, "b"],
                                            rotation = ellipse.positions[i, "rotation"],
                                            gp = gpar(lty = 0, fill = fill[ellipse.positions[i, "fill.mapping"]], alpha = alpha[ellipse.positions[i,"fill.mapping"]])))
  }
  for (i in 1:4) {
    grob.list <- gList(grob.list, ellipse(x = ellipse.positions[i, 
                                                                "x"], y = ellipse.positions[i, "y"], a = ellipse.positions[i, 
                                                                                                                           "a"], b = ellipse.positions[i, "b"], rotation = ellipse.positions[i, 
                                                                                                                                                                                             "rotation"], gp = gpar(lwd = lwd[ellipse.positions[i, 
                                                                                                                                                                                                                                                "line.mapping"]], lty = lty[ellipse.positions[i, 
                                                                                                                                                                                                                                                                                              "line.mapping"]], col = col[ellipse.positions[i, 
                                                                                                                                                                                                                                                                                                                                            "line.mapping"]], fill = "transparent")))
  }
  label.matrix <- matrix(nrow = 15, ncol = 3)
  colnames(label.matrix) <- c("label", "x", "y")
  label.matrix[1, ] <- c(a1, 0.35, 0.77)
  label.matrix[2, ] <- c(a2, 0.5, 0.69)
  label.matrix[3, ] <- c(a3, 0.65, 0.77)
  label.matrix[4, ] <- c(a4, 0.31, 0.67)
  label.matrix[5, ] <- c(a5, 0.4, 0.58)
  label.matrix[6, ] <- c(a6, 0.5, 0.47)
  label.matrix[7, ] <- c(a7, 0.6, 0.58)
  label.matrix[8, ] <- c(a8, 0.69, 0.67)
  label.matrix[9, ] <- c(a9, 0.18, 0.58)
  label.matrix[10, ] <- c(a10, 0.32, 0.42)
  label.matrix[11, ] <- c(a11, 0.425, 0.38)
  label.matrix[12, ] <- c(a12, 0.575, 0.38)
  label.matrix[13, ] <- c(a13, 0.68, 0.42)
  label.matrix[14, ] <- c(a14, 0.82, 0.58)
  label.matrix[15, ] <- c(a15, 0.5, 0.28)
  label.matrix[,1] <- round(label.matrix[,1],2);
  for (i in 1:nrow(label.matrix)) {
    grob.list <- gList(grob.list, textGrob(label = label.matrix[i,"label"],
                                           x = label.matrix[i, "x"],
                                           y = label.matrix[i, "y"],
                                           gp = gpar(col = label.col[i], cex = cex[i],
                                                     fontface = fontface[i], fontfamily = fontfamily[i])))
  }
  cat.pos.x <- c(0.18, 0.82, 0.35, 0.65)
  cat.pos.y <- c(0.58, 0.58, 0.77, 0.77)
  for (i in 1:4) {
    this.cat.pos <- find.cat.pos(x = cat.pos.x[i], y = cat.pos.y[i],
                                 pos = cat.pos[i], dist = cat.dist[i])
    grob.list <- gList(grob.list, textGrob(label = category[i],
                                           x = this.cat.pos$x, y = this.cat.pos$y, just = cat.just[[i]],
                                           gp = gpar(col = cat.col[i], cex = cat.cex[i], fontface = cat.fontface[i],
                                                     fontfamily = cat.fontfamily[i])))
  }
  grob.list <- VennDiagram::adjust.venn(VennDiagram::rotate.venn.degrees(grob.list,
                                                                         rotation.degree, rotation.centre[1], rotation.centre[2]),
                                        ...)
  if (ind) {
    grid.draw(grob.list)
  }
  return(grob.list)
}
myEnv <- new.env(parent = parent.frame());

if(!require(VennDiagram)){
  stop("Error: Package 'VennDiagram' cannot be loaded");
}
if((class(sets) != "list") || (length(sets) < 1)){
  stop("Error: 'sets' must be a non-empty list of vectors");
}
if(length(sets) > 5){
  stop("Error: cannot create a diagram with more than 5 sets");
}
args.extra <- as.list(match.call(expand.dots=TRUE));
args.extra[[1]] <- NULL;
for(arg in c("sets","add")){
  args.extra[[arg]] <- NULL;
}
ns <- length(sets);
nc <- length(unique(unlist(sets)));
## work out appropriate function
venn.fname <-
  c("single","pairwise","triple","quad","quintuple")[ns];
venn.fname <- sprintf("draw.%s.venn", venn.fname);
venn.func <- match.fun(venn.fname);
venn.df <- data.frame(row.names = unique(unlist(sets)));
for(x in 1:ns){
  venn.df[[sprintf("s%d",x)]] <- (rownames(venn.df) %in% sets[[x]]);
}
combs <- unlist(sapply(1:ns,function(m){combn(ns,m, simplify=FALSE)}),
                recursive=FALSE);
counts <- sapply(combs,function(x){
  sum(apply(as.matrix(venn.df[,x]),1,all));
});
if(percent){
  counts <- ((counts / nc) * 100);
}
## create names for count list
names(counts) <- sub("^n(.)$","area\\1",
                     sprintf("n%s",sapply(combs,paste,collapse="")));
oldnc <- names(counts);
## fiddle with names for single/pairwise to fit functions
if(ns == 1){
  names(counts) <- sub("area1","area",names(counts));
}
if(ns == 2){
  names(counts) <- sub("n12","cross.area",names(counts));
}
## replace colours of 'hack' plot so that regions can be identified
args.extraCol <- args.extra;
args.extraCol[["label.col"]] <- 1:length(oldnc);
args.extraCol[["cat.col"]] <- "black";
arg.list <- append(as.list(counts), args.extra);
arg.list$ind = FALSE; ## use 'add' for drawing plot
arg.listCol <- append(as.list(counts), args.extraCol);
arg.listCol$ind = FALSE; ## use 'add' for drawing plot
## hack to identify labels associated with particular regions
venn.plot.col <- do.call(venn.fname, arg.listCol, envir=parent.frame());
## actual (requested) plot
venn.plot <- do.call(venn.fname, arg.list, envir=parent.frame());
if(add == FALSE){
  grid.newpage();
}
grid.draw(venn.plot);
ven.res <- NULL;
sapply(venn.plot.col,function(obj){if(("text" %in% class(obj)) && (obj$gp$col != "black")){
  ven.res <<- rbind(ven.res,
                    data.frame(x=as.numeric(obj$x) ,y=as.numeric(obj$y), label=as.character(obj$label),
                               areaID=obj$gp$col))}});
## region -> grouping lookups generated from looking at the venn diagram source code
if(ns == 1){
  ven.res$columns <- c("1")[ven.res$areaID];
}
if(ns == 2){
  ven.res$groups <- c("1","1,2","2")[ven.res$areaID];
}
if(ns == 3){
  ven.res$groups <- c("1","1,2","2","1,3","1,2,3","2,3","3")[ven.res$areaID];
}
if(ns == 4){
  ven.res$groups <- c("3","3,4","4","1,3","1,3,4","1,2,3,4",
                      "2,3,4","2,4","1","1,4","1,2,4",
                      "1,2,3","2,3","2","1,2")[ven.res$areaID];
}
if(ns == 5){
  ven.res$groups <- c("1","2","3","4","5",
                      "3,5","1,5","1,4","1,2","2,5",
                      "2,3","1,3","3,4","2,4","4,5",
                      "3,4,5","1,3,5","1,4,5","1,2,4","1,2,5",
                      "2,3,5","1,2,3","1,3,4","2,3,4","2,4,5",
                      "2,3,4,5","1,3,4,5","1,2,4,5","1,2,3,5","1,2,3,4",
                      "1,2,3,4,5");
}
ven.res;
}