test_rbind_Assays <- function() {
    ## unnamed -- map by position
    l1 <- list(matrix(1, 3, 4), matrix(2, 3, 4))
    l2 <- list(matrix(3, 3, 4), matrix(4, 3, 4))
    a1 <- Assays(l1)
    a2 <- Assays(l2)
    current <- rbind(a1, a2)

    target <- Map(rbind, l1, l2)
    checkTrue(is(current, "ShallowSimpleListAssays"))
    checkIdentical(as(target, "SimpleList"), as(current, "SimpleList"))

    ## named -- map by name
    l1 <- list(x=matrix(1, 3, 4), y=matrix(2, 3, 4))
    l2 <- list(y=matrix(4, 3, 4), x=matrix(3, 3, 4))
    a1 <- Assays(l1)
    a2 <- Assays(l2)
    current <- rbind(a1, a2)

    target <- Map(rbind, l1, l2[match(names(l2), names(l1))])
    checkTrue(is(current, "ShallowSimpleListAssays"))
    checkIdentical(as(target, "SimpleList"), as(current, "SimpleList"))
}
