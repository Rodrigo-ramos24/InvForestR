#' Función de suma con omisión de NA
#'
#' @param x Vector o lista para suma con omision de NA
#'
#' @return El resultado de la suma con omision de NA
#' @export
#'
#' @examples
#' sum_na(x = c(34, 20, NA, 2))

sum_na <- function(x){
  c(yi = sum(x, na.rm = TRUE))
}


#' Estimador de razón por formaciones forestales para el INFyS
#'
#' @param yi Suma de la variable de interés
#' @param ai Suma de la variable auxiliar
#' @param Estrato Vector de estratos para la estimación del INFyS
#' @param Conglomerado Lista de conglomerados
#' @param AreasEstratos Vector de dos entradas: nombre del estrato y su superficie en ha
#'
#' @return Estimadores de razón para los estratos correspondientes
#'
#' @import doBy
#' @import stats
#'
#' @export
#'
#' @examples
#' # Datos de ejemplo
#' yi <- c(100, 200, 150, NA, 300)          # Sumas de la variable de interés
#' ai <- c(50, 75, 60, NA, 100)             # Sumas de la variable auxiliar
#' Estrato <- c("Estrato1", "Estrato2", "Estrato1", "Estrato2", "Estrato3")
#' Conglomerado <- c(1, 1, 2, 2, 3)        # Identificadores de conglomerado
#' AreasEstratos <- data.frame(Estrato = c("Estrato1", "Estrato2", "Estrato3"),
#'                              AreaHa = c(10, 20, 15)) # Superficies en ha
#'
#' # Llamada a la función
#' resultado <- ERaz(yi, ai, Estrato, Conglomerado, AreasEstratos)
#' print(resultado)
#'
ERaz <- function (yi, ai, Estrato, Conglomerado, AreasEstratos)
{
  n = length(yi)
  contadorSitios = rep(1, n)
  EstCong <- paste(as.character(Estrato), "-", as.character(Conglomerado))
  tmp = data.frame(yi, ai, contadorSitios, EstCong)
  BaseT1cong <- summaryBy(yi + ai + contadorSitios ~ EstCong,
                          data = tmp, FUN = sum_na, keep.names = TRUE, var.names = c("yi",
                                                                                     "ai",
                                                                                     "NumSitios"))
  BaseT1cong$Estrato <- substr(x = BaseT1cong$EstCong, start = 1,
                               stop = as.integer(gregexpr("-", BaseT1cong$EstCong)) -
                                 2)
  BaseT1cong$yi2 <- (BaseT1cong$yi)^2
  BaseT1cong$yiai <- (BaseT1cong$yi * BaseT1cong$ai)
  BaseT1cong$ai2 <- (BaseT1cong$ai)^2
  BaseT1cong$contadorCong <- 1
  BaseEstrato <- 0
  BaseEstrato <- summaryBy(NumSitios + contadorCong + yi +
                             ai + yi2 + yiai + ai2 ~ Estrato, data = BaseT1cong, FUN = sum_na,
                           keep.names = TRUE, var.names = c("NumSitios", "NumCong",
                                                            "yi", "ai", "yi2", "yiai", "ai2"))
  BaseEstrato <- merge(BaseEstrato, AreasEstratos, by.x = "Estrato",
                       by.y = "Estrato", all = FALSE)
  BaseEstrato$ERaz <- BaseEstrato$yi/BaseEstrato$ai
  BaseEstrato$Prom_ai <- (BaseEstrato$ai/BaseEstrato$NumCong)
  BaseEstrato$f <- 0
  BaseEstrato$f <- BaseEstrato$NumCong/BaseEstrato$AreaHa
  BaseEstrato$Var <- ((1 - BaseEstrato$f)/(BaseEstrato$NumCong *
                                             (BaseEstrato$NumCong - 1) * BaseEstrato$Prom_ai^2)) *
    (BaseEstrato$yi2 - 2 * BaseEstrato$ERaz * BaseEstrato$yiai +
       BaseEstrato$ai2 * (BaseEstrato$ERaz)^2)
  BaseEstrato$SdERaz <- sqrt(((1 - BaseEstrato$f)/(BaseEstrato$NumCong *
                                                     (BaseEstrato$NumCong - 1) * BaseEstrato$Prom_ai^2)) *
                               (BaseEstrato$yi2 - 2 * BaseEstrato$ERaz * BaseEstrato$yiai +
                                  BaseEstrato$ai2 * (BaseEstrato$ERaz)^2))
  BaseEstrato$Li <- BaseEstrato$ERaz - 2*BaseEstrato$SdERaz
  BaseEstrato$Ls <- BaseEstrato$ERaz + 2*BaseEstrato$SdERaz
  BaseEstrato$Error_muestreo <- (BaseEstrato$SdERaz/BaseEstrato$ERaz)
  BaseEstrato$p95 <- ((BaseEstrato$SdERaz*qt(0.975, n-1))/BaseEstrato$ERaz)
  BaseEstrato$p99 <- ((BaseEstrato$SdERaz*qt(0.995, n-1))/BaseEstrato$ERaz)
  as.data.frame(BaseEstrato)
}
