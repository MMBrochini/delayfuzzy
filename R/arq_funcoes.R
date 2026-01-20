#' Calcula o atraso de notificação em semanas
#'
#' @param dados Base de dados com datas
#' @return Base com coluna de atraso
#' @export
calcular_atraso <- function(dados) {
  dados$atrasoSemanas <- floor(dados$diferenca_3 / 7)
  return(dados)
}

#' Cria matriz de atrasos por semana epidemiológica
#'
#' @param dados Base com atraso calculado
#' @return Matriz de atraso
#' @export
matriz_atraso <- function(dados) {
  x <- factor(dados$SEM_NOT, levels = 1:53)
  y <- factor(dados$atrasoSemanas, levels = 0:52)
  as.matrix(table(x, y))
}

#' Gera matriz cumulativa de atrasos
#'
#' @param matriz Matriz de atraso
#' @return Matriz cumulativa
#' @export
matriz_cumulativa <- function(matriz) {
  sum_matriz <- matrix(0, nrow = nrow(matriz), ncol = ncol(matriz))
  sum_matriz[,1] <- matriz[,1]
  for (n in 2:ncol(matriz)) {
    sum_matriz[,n] <- sum_matriz[,n-1] + matriz[,n]
  }
  sum_matriz
}

#' Calcula o estimador Fj
#'
#' @param C Matriz cumulativa
#' @param j Semana de atraso
#' @param I Total de semanas
#' @return Valor de Fj
#' @export
calc_Fj <- function(C, j, I = 52) {
  indices_i <- 1:(I - j)
  sum(C[indices_i, j+1]) / sum(C[indices_i, j])
}

#' Gera números fuzzy triangulares
#'
#' @param Fjs Vetor de Fj
#' @param ljs Vetor de larguras
#' @return Lista de números fuzzy
#' @export
fuzzy_Fj <- function(Fjs, ljs) {
  left <- Fjs - ljs
  right <- Fjs + ljs
  sapply(seq_along(Fjs), function(i) {
    FuzzyNumbers::TriangularFuzzyNumber(left[i], Fjs[i], right[i])
  })
}
