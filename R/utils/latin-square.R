
# adapted from the JavaScript code at https://cs.uwaterloo.ca/~dmasson/tools/latin_square/
# which itself adapts Bradley (1958) algorithm
get_balanced_latin_square_order <- function (n_blocks, order_num) {
  
  # Bradley's method only holds for even numbers of conditions
  # and of course there are only n_blocks possible order numbers
  stopifnot(n_blocks %% 2 == 0, order_num <= n_blocks)
  
  # populate the first order by assigning conditions to odd positions 1, 3, 5, etc.
  # then backwards through the even numbered positions
  order <- rep(NA, n_blocks)
  order[seq(1, n_blocks, 2)] <- 1:(n_blocks/2)
  order[rev(seq(2, n_blocks, 2))] <- (n_blocks/2+1):n_blocks
  
  if (order_num == 1) return (order)
  
  # otherwise, add to each column position
  order <- order + (order_num - 1)
  # then circular increment
  # where if adding 1 to the highest condition num overflows it
  # replace it with the remainder to cycle back to
  order <- if_else(order > n_blocks,
                   order %% n_blocks,
                   order)
  
  return (order)
}

