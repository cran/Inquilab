#' Calculate First Order Kinetics Parameters
#'
#' This function calculates the rate constant, half-life, and provides a summary
#' of the first-order dissipation kinetics of pesticides, including the intercept,
#' R^2 value, and statistical measures of the fitted model.
#'
#' @param t Numeric vector, time points.
#' @param c Numeric vector, concentrations corresponding to each time point.
#'
#' @return A list containing the rate constant, half-life, and summary statistics
#' of the first order kinetics model.
#' @examples
#' t <- c(0, 5, 10, 15, 20, 25)
#' c <- c(100, 80, 60, 40, 20, 10)
#' first_order_kinetics(t, c)
#' @export
first_order_kinetics <- function(t, c) {
  # Ensure the data is in the expected format
  if(any(c <= 0)) stop("Concentrations must be greater than 0 to take the log")

  # Fit a first-order kinetic model
  kinetic_model <- lm(log(c) ~ t)

  # Extract the slope (negative of the rate constant) and calculate half-life
  rate_constant <- -coef(kinetic_model)[2]
  half_life <- log(2) / rate_constant

  # Summary of the model
  model_summary <- summary(kinetic_model)

  # Return the rate constant, half-life, and summary statistics
  list(rate_constant = rate_constant, half_life = half_life, summary = model_summary)
}



#' Calculate Second Order Kinetics Parameters
#'
#' This function calculates the rate constant and half-life based on second-order dissipation kinetics of pesticides,
#' and provides a summary of the kinetic model including intercept, R-squared value, and other statistical measures.
#'
#' @param t Numeric vector, time points.
#' @param c Numeric vector, concentrations corresponding to each time point.
#'
#' @return A list containing the rate constant(k), half-life  (t_half) for the second order kinetics, and a summary of the kinetic model.
#' @examples
#' t <- c(0, 5, 10, 15, 20, 25)
#' c <- c(100, 80, 60, 40, 20, 10)
#' result <- second_order_kinetics(t, c)
#' print(result$summary)
#' @export
second_order_kinetics <- function(t, c) {
  # Ensure the data is in the expected format
  if (any(c <= 0)) stop("Concentrations must be greater than 0 for calculation.")
  if (length(t) != length(c)) stop("Time and concentration vectors must be of equal length.")

  # Transform concentration for second order kinetics
  transformed_c <- 1 / c

  # Fit a linear model to the transformed data
  kinetic_model <- lm(transformed_c ~ t)

  # Calculate the rate constant and half-life for the second-order reaction
  rate_constant <- coef(kinetic_model)[2] # Corrected to kinetic_model
  half_life <- 1 / (rate_constant * min(c)) # Assuming initial concentration is the minimum value

  # Summary of the model
  model_summary <- summary(kinetic_model)

  # Return the rate constant, half-life, and summary statistics
  list(rate_constant = rate_constant, half_life = half_life, summary = model_summary)
}



#' Plot for First Order Kinetics
#'
#' This function plots the actual and predicted concentrations based on first-order kinetics.
#'
#' @param t Numeric vector, time points.
#' @param c Numeric vector, concentrations corresponding to each time point.
#' @param kinetic_model Model object, result of lm function fitting log(c) ~ t.
#'
#' @return Plot of actual vs. predicted concentrations for first order kinetics.
#' @examples
#' t <- c(0, 5, 10, 15, 20, 25)
#' c <- c(100, 80, 60, 40, 20, 10)
#' model <- lm(log(c) ~ t)
#' plot_first_order_kinetics(t, c, model)
#' @export
plot_first_order_kinetics <- function(t, c, kinetic_model) {
  plot(t, c, type = "o", xlab = "Time", ylab = "log (Concentration)",
       col = "blue", main = "First Order Kinetics")

  # Add the fitted line (predicted values) to the plot
  lines(t, exp(predict(kinetic_model)), col = "red", type = "o", lty = 4, lwd = 1)

  # Customize time axis breaks and labels (adjust as needed)
  axis.breaks <- c(0, 3, 7, 10, 15, 30, 45, 60, 90)
  axis(side = 1, at = axis.breaks, labels = axis.breaks)

}



#' Plot for Second Order Kinetics
#'
#' This function plots the actual and transformed (1/c) concentrations based on second-order kinetics.
#'
#' @param t Numeric vector, time points.
#' @param c Numeric vector, concentrations corresponding to each time point.
#' @param linear_model Model object, result of lm function fitting 1/c ~ t.
#'
#' @return Plot of actual vs. transformed concentrations for second order kinetics.
#' @examples
#' t <- c(0, 5, 10, 15, 20, 25)
#' c <- c(100, 80, 60, 40, 20, 10)
#' model <- lm(1/c ~ t)
#' plot_second_order_kinetics(t, c, model)
#' @export
plot_second_order_kinetics <- function(t, c, kinetic_model) {
  plot(t, 1 / c, type = "o", xlab = "Time", ylab = "1/(Concentration)",
       col = "blue", main = "Second Order Kinetics")

  # Add the fitted line (predicted values) to the plot
  lines(t, predict(kinetic_model), col = "red", type = "o", lty = 4, lwd = 1)

  # Customize time axis breaks and labels (adjust as needed)
  axis.breaks <- c(0, 3, 7, 10, 15, 30, 45, 60, 90)
  axis(side = 1, at = axis.breaks, labels = axis.breaks)

}



















