\name{Assays-class}
\docType{class}

\alias{class:Assays}
\alias{Assays-class}
\alias{Assays}
\alias{updateObject,Assays-method}
\alias{length,Assays-method}
\alias{names,Assays-method}
\alias{names<-,Assays-method}
\alias{getListElement,Assays-method}
\alias{setListElement,Assays-method}
\alias{dim,Assays-method}
\alias{[,Assays,ANY-method}
\alias{[,Assays,ANY,ANY,ANY-method}
\alias{[<-,Assays,ANY,ANY,ANY-method}
\alias{rbind,Assays-method}
\alias{cbind,Assays-method}
\alias{arbind,Matrix-method}
\alias{acbind,Matrix-method}

\alias{class:SimpleAssays}
\alias{SimpleAssays-class}
\alias{SimpleAssays}
\alias{coerce,SimpleList,SimpleAssays-method}
\alias{coerce,SimpleAssays,SimpleList-method}

\alias{class:ShallowData}
\alias{ShallowData-class}
\alias{ShallowData}

\alias{class:ShallowSimpleListAssays}
\alias{ShallowSimpleListAssays-class}
\alias{ShallowSimpleListAssays}
\alias{coerce,SimpleList,ShallowSimpleListAssays-method}
\alias{coerce,ShallowSimpleListAssays,SimpleList-method}

\alias{class:AssaysInEnv}
\alias{AssaysInEnv-class}
\alias{AssaysInEnv}
\alias{length,AssaysInEnv-method}
\alias{names,AssaysInEnv-method}
\alias{names<-,AssaysInEnv-method}
\alias{getListElement,AssaysInEnv-method}
\alias{setListElement,AssaysInEnv-method}
\alias{coerce,SimpleList,AssaysInEnv-method}
\alias{coerce,AssaysInEnv,SimpleList-method}

\title{Assays objects}

\description{
  The Assays virtual class and its methods provide a formal abstraction
  of the assays slot of \link{SummarizedExperiment} objects.

  SimpleAssays and ShallowSimpleListAssays are concrete subclasses of
  Assays with the former being currently the default implementation of Assays
  objects. Other implementations (e.g. disk-based) could easily be added.

  Note that these classes are not meant to be used directly by the end user
  and the material in this man page is aimed at package developers.
}

\details{
  Assays objects have a list-like semantics with elements having matrix- or
  array-like semantics (e.g., \code{dim}, \code{dimnames}).

  The Assays API consists of:
  \itemize{
    \item (a) The \code{Assays()} constructor function.
    \item (b) Lossless back and forth coercion from/to
              \link[S4Vectors]{SimpleList}. The coercion method from
              \link[S4Vectors]{SimpleList} doesn't need (and should not)
              validate the returned object.
    \item (c) \code{length}, \code{names}, \code{`names<-`},
              \code{getListElement}, \code{setListElement},
              \code{dim}, \code{[}, \code{`[<-`}, \code{rbind}, \code{cbind}.
  }
  An Assays concrete subclass needs to implement (b) (required) plus,
  optionally any of the methods in (c).

  IMPORTANT:
  \enumerate{
    \item Nobody in the Assays hierarchy is allowed to inherit from
          \link[S4Vectors]{SimpleList} because of the conflicting semantic
          of \code{[}.

    \item Methods that return a modified Assays object (a.k.a.
          endomorphisms), that is, \code{[} as well as replacement methods
          \code{names<-}, \code{setListElement}, and \code{[<-}, must respect
          the \emph{copy-on-change contract}.
          With objects that don't make use of references internally, the
          developer doesn't need to take any special action for that because
          it's automatically taken care of by R itself. However, for objects
          that do make use of references internally (e.g. environments,
          external pointers, pointer to a file on disk, etc...), the developer
          needs to be careful to implement endomorphisms with copy-on-change
          semantics. This can easily be achieved (and is what the default
          methods for Assays objects do) by performaing a full (deep) copy of
          the object before modifying it instead of trying to modify it
          in-place. However note that this full (deep) copy can be very
          expensive and is actually not necessary in order to achieve
          copy-on-change semantics: it's enough (and often preferrable for
          performance reasons) to copy only the parts of the object that need
          to be modified.
  }

  Assays has currently 3 implementations which are formalized by concrete
  subclasses SimpleAssays, ShallowSimpleListAssays, and AssaysInEnv.
  SimpleAssays is the default (prior to SummarizedExperiment 1.15.4,
  ShallowSimpleListAssays was the default). AssaysInEnv is a \emph{broken}
  alternative to ShallowSimpleListAssays that does NOT respect the
  \emph{copy-on-change contract}. It is only provided for illustration
  purposes (see source file Assays-class.R for the details).

  A little more detail about ShallowSimpleListAssays: a small reference
  class hierarchy (not exported from the \pkg{GenomicRanges} name space)
  defines a reference class ShallowData with a single field \code{data}
  of type \code{ANY}, and a derived class ShallowSimpleListAssays
  that specializes the type of \code{data} as \link[S4Vectors]{SimpleList},
  and \code{contains=c("ShallowData", "Assays")}. The assays slot of a
  \link{SummarizedExperiment} object contains an instance of
  ShallowSimpleListAssays.
}

\author{Martin Morgan and Hervé Pagès}

\seealso{
  \itemize{
    \item \link{SummarizedExperiment} objects.

    \item \link[S4Vectors]{SimpleList} objects in the \pkg{S4Vectors} package.
  }
}

\examples{
## ---------------------------------------------------------------------
## DIRECT MANIPULATION OF Assays OBJECTS
## ---------------------------------------------------------------------
m1 <- matrix(runif(24), ncol=3)
m2 <- matrix(runif(24), ncol=3)
a <- Assays(SimpleList(m1, m2))
a

as(a, "SimpleList")

length(a)
getListElement(a, 2)
dim(a)

b <- a[-4, 2]
b
length(b)
getListElement(b, 2)
dim(b)

names(a)
names(a) <- c("a1", "a2")
names(a)
getListElement(a, "a2")

rbind(a, a)
cbind(a, a)

## ---------------------------------------------------------------------
## COPY-ON-CHANGE CONTRACT
## ---------------------------------------------------------------------

## ShallowSimpleListAssays objects have copy-on-change semantics but not
## AssaysInEnv objects. For example:
ssla <- as(SimpleList(m1, m2), "ShallowSimpleListAssays")
aie <- as(SimpleList(m1, m2), "AssaysInEnv")

## No names on 'ssla' and 'aie':
names(ssla)
names(aie)

ssla2 <- ssla
aie2 <- aie
names(ssla2) <- names(aie2) <- c("A1", "A2")

names(ssla)  # still NULL (as expected)

names(aie)   # changed! (because the names<-,AssaysInEnv method is not
             # implemented in a way that respects the copy-on-change
             # contract)
}
