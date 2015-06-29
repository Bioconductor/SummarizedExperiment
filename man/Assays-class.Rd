\name{Assays-class}
\docType{class}

\alias{class:Assays}
\alias{Assays-class}
\alias{Assays}
\alias{length,Assays-method}
\alias{names,Assays-method}
\alias{names<-,Assays-method}
\alias{[[,Assays,ANY,ANY-method}
\alias{[[<-,Assays,ANY,ANY-method}
\alias{dim,Assays-method}
\alias{[,Assays,ANY-method}
\alias{[<-,Assays,ANY,ANY,ANY-method}
\alias{rbind,Assays-method}
\alias{cbind,Assays-method}

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
\alias{[[,AssaysInEnv,ANY,ANY-method}
\alias{coerce,SimpleList,AssaysInEnv-method}
\alias{coerce,AssaysInEnv,SimpleList-method}

\title{Assays objects}

\description{
  The Assays virtual class and its implementations (ShallowSimpleListAssays
  and AssaysInEnv) are used for storing the assays of a
  \link{SummarizedExperiment0} object. These classes are not meant to be
  used directly by the end-user.
  This man page contains material meant for package developers.
}

\details{
  Assays objects have a list-like semantics with elements having matrix- or
  array-like semantics (e.g., \code{dim}, \code{dimnames}).

  The Assays API consists of:
  \itemize{
    \item (a) the \code{Assays()} constructor function
    \item (b) lossless back and forth coercion from/to \link{SimpleList}
    \item (c) \code{length}, \code{names}, \code{`names<-`},
              \code{[[}, \code{`[[<-`},
              \code{dim}, \code{[}, \code{`[<-`}, \code{rbind}, \code{cbind}.
  }
  An Assays croncrete subclass needs to implement (b) (required) plus,
  optionally any of the methods in (c).

  IMPORTANT: Replacement methods (\code{names<-}, \code{[[<-}, \code{[<-})
  \emph{must} copy the object instead of trying to modify it in-place.
  The default methods for Assays objects do that. This ensures that Assays
  objects have a copy-on-change semantics. Note that this is particularly
  important for croncrete subclasses like ShallowSimpleListAssays that are
  based on objects with reference semantics.

  Assays has currently 2 implemenations which are formalized by concrete
  subclasses ShallowSimpleListAssays and AssaysInEnv, the former being the
  default.

  About ShallowSimpleListAssays: In a little more detail, a small
  reference class hierarchy (not exported from the \pkg{GenomicRanges}
  name space) defines a reference class ShallowData with a single field
  \code{data} of type \code{ANY}, and a derived class ShallowSimpleListAssays
  that specializes the type of \code{data} as \link{SimpleList}, and
  \code{contains=c("ShallowData", "Assays")}. The assays slot of a
  \link{SummarizedExperiment0} object contains an instance of
  ShallowSimpleListAssays.
}

\author{Martin Morgan, \url{mtmorgan@fhcrc.org}}

\seealso{
  \itemize{
    \item \link{SummarizedExperiment0} objects.

    \item \link[S4Vectors]{SimpleList} objects in the \pkg{S4Vectors} package.
  }
}

\examples{
m1 <- matrix(runif(24), ncol=3)
m2 <- matrix(runif(24), ncol=3)
a <- Assays(SimpleList(m1, m2))
a

length(a)
a[[2]]
dim(a)

b <- a[-4, 2]
b
length(b)
b[[2]]
dim(b)

names(a)
names(a) <- c("a1", "a2")
names(a)
a[["a2"]]

rbind(a, a)
cbind(a, a)
}
