
  # This code generates a color palette using the Leiden University colors
  # source: https://drsimonj.svbtle.com/creating-corporate-colour-palettes-for-ggplot2
  
  # Create vector with 8 Leiden university / LACDR colors
  lei_colors <- c(
    `blue`       = "#001158",
    `orange`     = "#FF9933", 
    `red`        = "#be1908",
    `lightgreen` = "#aaad00",
    `brightgreen`= "#06bd09",
    `darkgreen`  = "#2c712d",
    `turquoise`  = "#34a3a9",
    `lightblue`  = "#5cb1eb",
    `brightblue` = "#0536fa",
    `violet`     = "#b02079" )
  
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
    `nine`  = lei_cols("lightblue", "violet", "brightgreen",  "brightblue", "red", "lightgreen", "blue", "orange", "darkgreen"))
  
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
  
