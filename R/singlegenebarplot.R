# --------------
# Date:  2019-11-18 13:11:32
# Author:JijunYu
# Email: jijunyu140@gmail.com
# --------------
#' About project:This is for barplot
#' @param inputdata This is a data frame
#' @param x x varibles for drawing (cohort name,character)
#' @param y y varibles for drawing (values,character)
#' @param fill fill colors (could be color name)
#' @param color  ourline colors
#' @param sort.val 	a string specifying whether the value should be sorted. Allowed values are "none" (no sorting), "asc" (for ascending) or "desc" (for descending).
#' @param sort.by.groups 	logical value. If TRUE the data are sorted by groups. Used only when sort.val != "none".
#' @param  x.text.angle Rotate vertically x axis texts
#' @param ylab 	character vector specifying y axis labels. Use ylab = FALSE to hide ylab.
#' @param legend.title character vector for legend title
#' @param rotate rotate angle
#' @param ggtheme function, ggplot2 theme name. Default value is theme_pubr(). Allowed values include ggplot2 official themes: theme_gray(), theme_bw(), theme_minimal(), theme_classic(), theme_void(), ....
#' @export
singlegenebarplot <- function(inputdata, x,y,fill,sort.by.groups,ylab,xlab,legend.title,title){
  require("ggpubr")
  p <- ggbarplot(inputdata, x = x, y = y,
            fill = fill,           # change fill color by mpg_level
            color = "white",            # Set bar border colors to white
            palette = "jco",            # jco journal color palett. see ?ggpar
            sort.val = "desc",          # Sort the value in descending order
            sort.by.groups = sort.by.groups,     # Don't sort inside each group
            x.text.angle = 90,
            xlab = xlab,          # Rotate vertically x axis texts
            ylab = ylab,
            legend.title = legend.title,
            rotate = TRUE,
            ggtheme = theme_bw(),
            title = title
  ) +
    theme(legend.title = element_text(face = "bold")) +
    theme(title = element_text(face = "bold",size = 20,hjust = 0.5)) +
    theme(axis.text = element_text(face = "bold")) +
    theme(legend.position = "none") +
    theme(axis.text.y = element_text(face = "bold",size=15))
  plot(p)
}

