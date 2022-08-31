
# This code generates a color palette using the Leiden University colors
# source: https://drsimonj.svbtle.com/creating-corporate-colour-palettes-for-ggplot2

# Create vector with 8 Leiden university / LACDR colors
lei_colors <- c(
  `blue`        = "#001158",
  `orange`      = "#FF9933",
  `red`         = "#be1908",
  `lightgreen`  = "#aaad00",
  `brightgreen` = "#06bd09",
  `darkgreen`   = "#2c712d",
  `turquoise`   = "#34a3a9",
  `lightblue`   = "#5cb1eb",
  `brightblue`  = "#0536fa",
  `violet`      = "#b02079",
  `black`       = "#000000",
  `soft orange` = "#f46e32",
  `soft blue`   = "#8592BC")

## Function to extract the hex codes from this vector by name
#' Function to extract lei colors as hex codes
#' @param ... Character names of lei_colors 
lei_cols <- function (...){
  cols <- c(...)
  
  if (is.null(cols))
    return (lei_colors)
  
  lei_colors[cols]
}

lei_palettes <- list(
  `main`  = lei_cols("blue", "orange"),
  `three` = lei_cols("blue", "orange", "darkgreen"),
  `cool`  = lei_cols("blue", "lightblue", "turquoise", "lightgreen", "darkgreen"),
  `hot`   = lei_cols("violet", "red", "orange"),
  `mixed` = lei_cols("blue", "lightblue", "turquoise", "darkgreen", "lightgreen", "orange", "red", "violet"),
  `two`   = lei_cols("red", "violet"), 
  `five`  = lei_cols("blue", "lightblue", "orange", "red", "darkgreen"),
  `six`  = lei_cols("blue", "lightblue", "orange", "red", "darkgreen", "black"),
  `nine`  = lei_cols("lightblue", "violet", "brightgreen",  "brightblue", "red", "lightgreen", "blue", "orange", "darkgreen"),
  `gradient` = lei_cols("blue", "soft blue", "soft orange"),
  `gradient2` = lei_cols("blue", "soft blue", "orange", "soft orange"))

#' Return function to interpolate a lei color palette
#' @param palette Character name of palette in lei_palettes
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments to pass to colorRampPalette(), such as an alpha value
lei_pal <- function(palette = "main", reverse = FALSE, ...) {
  pal <- lei_palettes[[palette]]
  
  if (reverse) pal <- rev(pal)
  
  colorRampPalette(pal, ...)
}

# Now return a function for any palette, for example `cool`:
lei_pal("cool")
# The returned function will interpolate the palette colors for a certain number of levels, making it possible to create shades between our original colors. To demonstrate, we can interpolate the "cool" palette to a length of 10:
lei_pal("cool")(10)
# This is what we need to create custom ggplot2 scales

## Create custom color and fill scales for ggplot2 by creating one function for color and one for fill. 
#' Color scale constructor for lei colors
#' @param palette Character name of palette in lei_palettes
#' @param discrete Boolean indicating whether color aesthetic is discrete or not
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments passed to discrete_scale() or
#'            scale_color_gradientn(), used respectively when discrete is TRUE or FALSE
#'
scale_color_lei <- function(palette = "main", discrete = TRUE, reverse = FALSE, ...) {
  pal <- lei_pal(palette = palette, reverse = reverse)
  
  if (discrete) {
    discrete_scale("colour", paste0("lei_", palette), palette = pal, ...)
  } else {
    scale_color_gradientn(colours = pal(256), ...)
  }
}

#' Fill scale constructor for lei colors
#'
#' @param palette Character name of palette in lei_palettes
#' @param discrete Boolean indicating whether color aesthetic is discrete or not
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments passed to discrete_scale() or
#'            scale_fill_gradientn(), used respectively when discrete is TRUE or FALSE
#'
scale_fill_lei <- function(palette = "main", discrete = TRUE, reverse = FALSE, ...) {
  pal <- lei_pal(palette = palette, reverse = reverse)
  
  if (discrete) {
    discrete_scale("fill", paste0("lei_", palette), palette = pal, ...)
  } else {
    scale_fill_gradientn(colours = pal(256), ...)
  }
}

##Function for scaling x axis from 0:4 and 30
scale_x_time_squish <- function(ticks = 0:30, showing = c(0:4,  30), ...) {
  ticks[!ticks %in% showing] <- ""
  
  squish_trans <- function(from, to, factor) {
    trans <- function(x) {
      if (any(is.na(x))) return(x)
      # get indices for the relevant regions
      isq <- x > from & x < to
      ito <- x >= to
      # apply transformation
      x[isq] <- from + (x[isq] - from)/factor
      x[ito] <- from + (to - from)/factor + (x[ito] - to)
      return(x)
    }
    inv <- function(x) {
      if (any(is.na(x))) return(x)
      # get indices for the relevant regions
      isq <- x > from & x < from + (to - from)/factor
      ito <- x >= from + (to - from)/factor
      # apply transformation
      x[isq] <- from + (x[isq] - from) * factor
      x[ito] <- to + (x[ito] - (from + (to - from)/factor))
      return(x)
    }
    # return the transformation
    return(scales::trans_new("squished", trans, inv))
  }
  
  scale_x_continuous(trans = squish_trans(4.5, 29.5, 7), breaks = 0:30, labels = as.character(ticks), ...)
}




