#Function to get the layer sizes
get.layer.sizes <- function(input.data, output.data, hidden.neurons, train=TRUE) {
  n.X <- dim(input.data)[1]
  n.hidden <- hidden.neurons
  n.y <- dim(output.data)[1]   
  
  sizes <- list("n.X" = n.X,
                "n.hidden" = n.hidden,
                "n.y" = n.y)
  
  return(sizes)
}

#Function to initialize neural network parameters
initialize.parameters <- function(input.data, layer.sizes){
  m <- dim(data.matrix(input.data))[2]
  
  n.X <- layer.sizes$n.X
  n.hidden <- layer.sizes$n.hidden
  n.y <- layer.sizes$n.y
  
  #Initialize weight matrices and bias vectors with small random values
  weights.1 <- matrix(runif(n.hidden * n.X), nrow = n.hidden, ncol = n.X, byrow = TRUE) * 0.01
  biases.1 <- matrix(rep(0, n.hidden), nrow = n.hidden)
  weights.2 <- matrix(runif(n.y * n.hidden), nrow = n.y, ncol = n.hidden, byrow = TRUE) * 0.01
  biases.2 <- matrix(rep(0, n.y), nrow = n.y)
  
  parameters <- list("weights.1" = weights.1,
                     "biases.1" = biases.1, 
                     "weights.2" = weights.2,
                     "biases.2" = biases.2)
  
  return (parameters)
}

#Function to perform forward propagation through the neural network
forward.propagation <- function(input.data, parameters, layer.sizes){
  m <- dim(input.data)[2]
  n.hidden <- layer.sizes$n.hidden
  n.y <- layer.sizes$n.y
  
  weights.1 <- parameters$weights.1
  biases.1 <- parameters$biases.1
  weights.2 <- parameters$weights.2
  biases.2 <- parameters$biases.2
  
  #Calculate intermediate values in the network
  z.1 <- weights.1 %*% input.data + matrix(rep(biases.1, m), nrow = n.hidden)
  a.1 <- sigmoid(z.1)
  z.2 <- weights.2 %*% a.1 + matrix(rep(biases.2, m), nrow = n.y)
  a.2 <- sigmoid(z.2)
  
  cache <- list("z.1" = z.1,
                "a.1" = a.1, 
                "z.2" = z.2,
                "a.2" = a.2)
  
  return (cache)
}

#Function to compute the cost (error) of the neural network
compute.cost <- function(input.data, output.data, cache) {
  m <- dim(input.data)[2]
  a.2 <- cache$a.2
  
  #Calculate the cost using the binary cross-entropy loss
  logprobs <- (log(a.2) * output.data) + (log(1-a.2) * (1-output.data))
  cost <- -sum(logprobs/m)
  
  return (cost)
}

#Function to perform backward propagation and compute gradients
backward.propagation <- function(input.data, output.data, cache, parameters, layer.sizes){
  m <- dim(input.data)[2]
  n.X <- layer.sizes$n.X
  n.hidden <- layer.sizes$n.hidden
  n.y <- layer.sizes$n.y
  
  a.2 <- cache$a.2
  a.1 <- cache$a.1
  weights.2 <- parameters$weights.2
  
  #Calculate gradients using backpropagation
  dz.2 <- a.2 - output.data
  dw.2 <- 1/m * (dz.2 %*% t(a.1)) 
  db.2 <- matrix(1/m * sum(dz.2), nrow = n.y)
  db.2.new <- matrix(rep(db.2, m), nrow = n.y)
  
  dz.1 <- (t(weights.2) %*% dz.2) * (1 - a.1^2)
  dw.1 <- 1/m * (dz.1 %*% t(input.data))
  db.1 <- matrix(1/m * sum(dz.1), nrow = n.hidden)
  db.1.new <- matrix(rep(db.1, m), nrow = n.hidden)
  
  gradients <- list("dw.1" = dw.1, 
                    "db.1" = db.1,
                    "dw.2" = dw.2,
                    "db.2" = db.2)
  
  return(gradients)
}

#Function to update neural network parameters based on gradients
update.parameters <- function(gradients, parameters, learning.rate){
  weights.1 <- parameters$weights.1
  biases.1 <- parameters$biases.1
  weights.2 <- parameters$weights.2
  biases.2 <- parameters$biases.2
  
  dw.1 <- gradients$dw.1
  db.1 <- gradients$db.1
  dw.2 <- gradients$dw.2
  db.2 <- gradients$db.2
  
  #Update parameters using gradient descent
  weights.1 <- weights.1 - learning.rate * dw.1
  biases.1 <- biases.1 - learning.rate * db.1
  weights.2 <- weights.2 - learning.rate * dw.2
  biases.2 <- biases.2 - learning.rate * db.2
  
  updated.parameters <- list("weights.1" = weights.1,
                             "biases.1" = biases.1,
                             "weights.2" = weights.2,
                             "biases.2" = biases.2)
  
  return (updated.parameters)
}

#Function to train the neural network
train <- function(input.data, output.data, epochs, hidden.neurons, learning.rate){
  layer.sizes <- get.layer.sizes(input.data, output.data, hidden.neurons)
  initial.parameters <- initialize.parameters(input.data, layer.sizes)
  
  cost.history <- c()
  
  for (i in 1:epochs) {
    cache <- forward.propagation(input.data, initial.parameters, layer.sizes)
    cost <- compute.cost(input.data, output.data, cache)
    gradients <- backward.propagation(input.data, output.data, cache, initial.parameters, layer.sizes)
    updated.parameters <- update.parameters(gradients, initial.parameters, learning.rate = learning.rate)
    initial.parameters <- updated.parameters
    
    cost.history <- c(cost.history, cost)
    
    if (i %% 10000 == 0) cat("Iteration", i, " | Cost: ", cost, "\n")
  }
  
  result <- list("updated.parameters" = updated.parameters,
                 "cost.history" = cost.history)
  
  return (result)
}

