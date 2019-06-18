test_bind_Assays <- function() {
    ## unnamed -- map by position
    l1 <- list(matrix(1, 3, 4), matrix(2, 3, 4))
    l2 <- list(matrix(3, 3, 4), matrix(4, 3, 4))
    a1 <- Assays(l1)
    a2 <- Assays(l2)

    target <- Map(rbind, l1, l2)
    current <- rbind(a1, a2)
    checkTrue(is(current, "SimpleAssays"))
    checkIdentical(as(target, "SimpleList"), as(current, "SimpleList"))

    target <- Map(cbind, l1, l2)
    current <- cbind(a1, a2)
    checkTrue(is(current, "SimpleAssays"))
    checkIdentical(as(target, "SimpleList"), as(current, "SimpleList"))

    ## named -- map by name
    l1 <- list(x=matrix(1, 3, 4), y=matrix(2, 3, 4))
    l2 <- list(y=matrix(4, 3, 4), x=matrix(3, 3, 4))
    a1 <- Assays(l1)
    a2 <- Assays(l2)

    target <- Map(rbind, l1, l2[match(names(l2), names(l1))])
    current <- rbind(a1, a2)
    checkTrue(is(current, "SimpleAssays"))
    checkIdentical(as(target, "SimpleList"), as(current, "SimpleList"))

    target <- Map(cbind, l1, l2[match(names(l2), names(l1))])
    current <- cbind(a1, a2)
    checkTrue(is(current, "SimpleAssays"))
    checkIdentical(as(target, "SimpleList"), as(current, "SimpleList"))
}

test_bind_higher_order_Assays <- function() {
    ## unnamed -- map by position
    l1 <- list(array(1, dim = c(3, 4, 5, 6)),
               array(2, dim = c(3, 4, 5, 6)))
    l2 <- list(array(4, dim = c(3, 4, 5, 6)),
               array(3, dim = c(3, 4, 5, 6)))
    a1 <- Assays(l1)
    a2 <- Assays(l2)

    target <- Map(arbind, l1, l2)
    current <- rbind(a1, a2)
    checkTrue(is(current, "SimpleAssays"))
    checkIdentical(as(target, "SimpleList"), as(current, "SimpleList"))

    target <- Map(acbind, l1, l2)
    current <- cbind(a1, a2)
    checkTrue(is(current, "SimpleAssays"))
    checkIdentical(as(target, "SimpleList"), as(current, "SimpleList"))

    ## named -- map by name
    l1 <- list(x = array(1, dim = c(3, 4, 5, 6)),
               y = array(2, dim = c(3, 4, 5, 6)))
    l2 <- list(y = array(4, dim = c(3, 4, 5, 6)),
               x = array(3, dim = c(3, 4, 5, 6)))
    a1 <- Assays(l1)
    a2 <- Assays(l2)

    target <- Map(arbind, l1, l2[match(names(l2), names(l1))])
    current <- rbind(a1, a2)
    checkTrue(is(current, "SimpleAssays"))
    checkIdentical(as(target, "SimpleList"), as(current, "SimpleList"))

    target <- Map(acbind, l1, l2[match(names(l2), names(l1))])
    current <- cbind(a1, a2)
    checkTrue(is(current, "SimpleAssays"))
    checkIdentical(as(target, "SimpleList"), as(current, "SimpleList"))
}

test_bind_error_on_incompatible_dimension_Assays <- function() {
    l1 <- list(x = array(1, dim = c(3, 4, 5, 6)),
               y = array(2, dim = c(3, 4, 5, 6)))
    l2 <- list(y = matrix(4, 3, 4), x = matrix(3, 3, 4))
    a1 <- Assays(l1)
    a2 <- Assays(l2)

    ## arbind
    checkException(rbind(a1, a2), silent = TRUE)
    checkException(arbind(l1[[1]], l2[[1]]), silent = TRUE)

    ## acbind
    checkException(cbind(a1, a2), silent = TRUE)
    checkException(acbind(l1[[1]], l2[[1]]), silent = TRUE)
}

