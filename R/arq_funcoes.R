#' Calcula o atraso de notificação em semanas
#'
#' @param dados Base de dados
#' @param data_notificacao Coluna da data de notificação
#' @param data_sintomas Coluna da data de início dos sintomas
#'
#' @return Base com coluna atrasoSemanas
#' @export
calcular_atraso <- function(dados, data_notificacao, data_sintomas) {

  diferenca_dias <- as.Date(dados[[data_notificacao]]) -
    as.Date(dados[[data_sintomas]])

  dados$atrasoSemanas <- floor(as.numeric(diferenca_dias) / 7)

  return(dados)
}


#' Cria matriz de atrasos por semana epidemiológica
#'
#' @param dados Base com atraso calculado
#' @return Matriz de atraso (frequências)
#' @export
matriz_atraso <- function(dados) {

# contagem total de registros
  n_total <- nrow(dados)

# identificar registros válidos
  validos <- !is.na(dados$SEM_NOT) &
    !is.na(dados$atrasoSemanas) &
    dados$SEM_NOT %in% 1:53 &
    dados$atrasoSemanas %in% 0:52

  n_validos <- sum(validos)
  n_excluidos <- n_total - n_validos

  if (n_excluidos > 0) {
    warning(
      paste0(
        "Foram excluídos ", n_excluidos, " registros na construção da matriz de atraso. ",
        "Apenas casos com semana epidemiológica válida (1–53) ",
        "e atraso de notificação entre 0 e 52 semanas foram considerados. ",
        "Registros com atrasos negativos, atrasos extremos ou datas inconsistentes ",
        "não são incluídos na matriz."
      ),
      call. = FALSE
    )
  }

  x <- factor(dados$SEM_NOT[validos], levels = 1:53)
  y <- factor(dados$atrasoSemanas[validos], levels = 0:52)

  table(x, y)
}


#' Gera matriz cumulativa de atrasos
#'
#' @param matriz Matriz de atraso
#' @return Matriz cumulativa
#' @export
matriz_cumulativa <- function(matriz) {
  sum_matriz <- matrix(0, nrow = nrow(matriz), ncol = ncol(matriz))
  sum_matriz[, 1] <- matriz[, 1]
  for (n in 2:ncol(matriz)) {
    sum_matriz[, n] <- sum_matriz[, n - 1] + matriz[, n]
  }
  sum_matriz
}

#' Estima o fator de correção Fj
#'
#' @param C Matriz cumulativa de atrasos.
#' @param j Índice do atraso (semana).
#'
#' @return Valor escalar do estimador Fj.
#' @export
estimar_Fj <- function(C, j) {

  I <- nrow(C)

  if (I - j - 1 < 1) {
    return(NA_real_)
  }

  indices_i <- 1:(I - j - 1)

  sum(C[indices_i, j], na.rm = TRUE) /
    sum(C[indices_i, j + 1], na.rm = TRUE)
}


#' Calcula a largura fuzzy lj associada ao estimador Fj
#'
#' @param C Matriz cumulativa.
#' @param X Matriz de atraso (frequências).
#' @param j Índice do atraso.
#'
#' @return Valor escalar lj.
#' @export
calc_lj <- function(C, X, j) {

  I <- nrow(C)

  if (I - j - 1 < 1) {
    return(NA_real_)
  }

  indices_i <- 1:(I - j - 1)

  sum(X[indices_i, j + 1], na.rm = TRUE) /
    sum(C[indices_i, j], na.rm = TRUE)
}



#' Multiplicação de números fuzzy triangulares
#'
#' @param A Vetor (centro, esquerda, direita).
#' @param B Vetor (centro, esquerda, direita).
#'
#' @return Vetor fuzzy triangular resultante.
#' @export
mult_fuzzy <- function(A, B) {

  a <- A[1]; la <- A[2]; ra <- A[3]
  b <- B[1]; lb <- B[2]; rb <- B[3]

  c(
    a * b,
    a * lb + b * la - la * lb,
    a * rb + b * ra + ra * rb
  )
}


#' Propaga incerteza fuzzy ao longo da matriz cumulativa
#'
#' @param C Matriz cumulativa.
#' @param Fjs Vetor de fatores Fj.
#' @param ljs Vetor de larguras fuzzy lj.
#'
#' @return Lista com matrizes centro, esquerda e direita.
#' @export
propagar_fuzzy <- function(C, Fjs, ljs) {

  I <- nrow(C)
  J <- ncol(C)

  Cc <- C
  left_C  <- matrix(0, nrow = I, ncol = J)
  right_C <- matrix(0, nrow = I, ncol = J)

  for (i in 2:I) {
    for (j in (I - i + 1):(I - 1)) {

      if (j + 1 > J || is.na(Fjs[j]) || is.na(ljs[j])) next

      Cij <- c(Cc[i, j], left_C[i, j], right_C[i, j])
      Fj  <- c(Fjs[j], ljs[j], ljs[j])

      Cij1 <- mult_fuzzy(Fj, Cij)

      Cc[i, j + 1]       <- Cij1[1]
      left_C[i, j + 1]   <- Cij1[2]
      right_C[i, j + 1]  <- Cij1[3]
    }
  }

  list(
    centro = Cc,
    esquerda = left_C,
    direita = right_C
  )
}



#' Correção de atraso por estimadores fuzzy
#'
#' @param dados Base de dados.
#' @param data_notificacao Coluna de notificação.
#' @param data_sintomas Coluna de início dos sintomas.
#'
#' @return Lista com matrizes intermediárias, estimadores e previsões fuzzy.
#' @export
estimacao_delay_fuzzy <- function(dados,
                                  data_notificacao,
                                  data_sintomas) {

  # 1. Calcular atraso em semanas
  dados <- calcular_atraso(dados, data_notificacao, data_sintomas)

  # 2. Matriz de atrasos observados
  X <- matriz_atraso(dados)

  # 3. Matriz cumulativa
  C <- matriz_cumulativa(X)

  # 4. Definição correta do domínio de j
  I <- nrow(C)
  j_validos <- 1:(I - 2)

  # 5. Estimadores Fj e larguras lj
  Fjs <- sapply(j_validos, estimar_Fj, C = C)
  ljs <- sapply(j_validos, calc_lj, C = C, X = X)

  # 6. Propagação fuzzy
  previsao <- propagar_fuzzy(C, Fjs, ljs)

  list(
    matriz_atraso     = X,
    matriz_cumulativa = C,
    Fj                = Fjs,
    lj                = ljs,
    previsao          = previsao
  )
}
