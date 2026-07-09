# Fortran Programming Paradigms — A Practical Guide for Computational Mathematics

> Written for: a doctoral student in computational mathematics (radio spectrum mapping), intermediate Fortran proficiency, using fortran-lang/stdlib and `fpm`.

---

## 1. Paradigm Quick Reference

| Paradigm | Core Idea | Fortran Support | Verdict |
|---|---|---|---|
| **Structured** | Procedures (subroutines/functions), clear control flow, no `GOTO` spaghetti | ✅ Native — the baseline of modern Fortran | Use everywhere |
| **Modular** | Code organized into `module` units with explicit `public`/`private` interfaces | ✅ First-class since Fortran 90 | **Primary** |
| **Array** | Operations on whole arrays: `C = matmul(A, B)`, `where`, array sections | ✅ Fortran's killer feature | **Primary** |
| **Object-Oriented** | Derived types with type-bound procedures, inheritance, polymorphism | ⚠️ Supported since F2003, verbose; no GC | **Use selectively** |
| **Generic** | Algorithms written once, work for `real32`, `real64`, `complex`, etc. | ⚠️ `template` in F2023 (new, limited compiler support) | Wait for compiler maturity |
| **Functional** | Pure functions, no side effects, immutable state | ⚠️ Partial — `pure`/`elemental` exist, but no closures or algebraic data types | Skip as a primary style |

---

## 2. Paradigm-to-Technique Mapping

| Technique | Natural Paradigm | Why |
|---|---|---|
| **Numerical linear algebra** | **Array** | Fortran's array intrinsincs map directly to mathematical notation |
| **Gaussian process regression** | **Array + Modular** | Kernels produce matrices → array ops; kernel variants → module organization |
| **Bayesian inference** | **Structured + Modular** | Samplers are algorithmic pipelines; modules by sampler/distribution family |
| **Optimization** | **OOP + Structured** | Different optimizers share a common interface — OOP earns its keep here |
| **Radio spectrum mapping** | **Array + Modular** | Heavy linear algebra on grids; module organization by method |

---

## 3. Recommendation: A Pragmatic Hybrid

### Core stack (70–80% of code): **Modular + Array**

```fortran
! ─── mod_kernels.f90 ───
module mod_kernels
  implicit none
  private
  public :: rbf_kernel, matern32_kernel

contains

  ! Elemental: operates element-wise on arrays of any rank
  elemental real function rbf_kernel(x, y, lengthscale) result(k)
    real, intent(in) :: x, y, lengthscale
    k = exp(-0.5 * ((x - y) / lengthscale)**2)
  end function

  ! Pure: no side effects, safe for parallelization
  pure real function matern32_kernel(x, y, lengthscale) result(k)
    real, intent(in) :: x, y, lengthscale
    real :: r
    r = sqrt(3.0) * abs(x - y) / lengthscale
    k = (1.0 + r) * exp(-r)
  end function
end module
```

```fortran
! ─── mod_gp.f90 ───
module mod_gp
  use mod_kernels, only: rbf_kernel
  implicit none
  private
  public :: gp_predict

contains

  ! Array programming: matrix-level operations, not element loops
  function gp_predict(X_train, y_train, X_test, lengthscale, noise) result(mu)
    real, intent(in) :: X_train(:,:), y_train(:), X_test(:,:)
    real, intent(in) :: lengthscale, noise
    real, allocatable :: mu(:), K(:,:), K_star(:,:), L(:,:), alpha(:)
    integer :: n

    n = size(X_train, 1)

    ! Build covariance matrix — array operations
    K = build_covariance_matrix(X_train, X_train, lengthscale)
    K = K + noise * eye(n)                       ! add noise to diagonal
    K_star = build_covariance_matrix(X_test, X_train, lengthscale)

    ! Cholesky + solve
    call potrf(K)                                  ! in-place Cholesky
    alpha = solve(K, y_train)                      ! forward/back substitution
    mu = matmul(K_star, alpha)                     ! predictive mean
  end function
end module
```

### Selective OOP: Only for interchangeable components with 3+ variants

- **Optimizers** (L-BFGS, Adam, conjugate gradient, trust-region)
- **Kernels** (RBF, Matérn, periodic, rational quadratic)
- **Likelihoods** (Gaussian, Student-t, Poisson)
- **Samplers** (Metropolis-Hastings, HMC, NUTS, slice sampling)
- **Preconditioners** (Jacobi, ILU, AMG)

---

## 4. The Paradigm Decision Tree

```
┌─ Is the math expressed as matrix/vector equations?
│  └─ YES → Array programming
│
├─ Do I need to swap between 3+ variants of the same thing?
│  └─ YES → OOP with abstract type
│
├─ Is it a standalone algorithm with clear steps?
│  └─ YES → Structured/modular
│
├─ Is it a utility that works element-by-element on scalars and arrays?
│  └─ YES → elemental function
│
└─ Default → Module + functions with array arguments
```

---

## 5. Plain Derived Types vs. OOP Derived Types

### Level 0 — Plain Derived Type (Just Data)

A C-style struct. The type is passive — just a data container.

```fortran
module mod_plain
  implicit none
  private
  public :: gp_params, fit_gp, predict_gp

  type :: gp_params
    real    :: lengthscale = 1.0
    real    :: noise       = 1e-6
    integer :: n_iter      = 100
  end type

contains

  subroutine fit_gp(self, X, y)
    type(gp_params), intent(inout) :: self
    real, intent(in) :: X(:,:), y(:)
  end subroutine

  function predict_gp(self, X_test) result(mu)
    type(gp_params), intent(in) :: self
    real, intent(in) :: X_test(:,:)
    real, allocatable :: mu(:)
  end function
end module
```

**Use when:** The type is a simple parameter bundle or data record (config structs, point coordinates, dataset headers). No behavioral variation.

---

### Level 1 — Type-Bound Procedures (Encapsulation, No Inheritance)

Procedures attached to the type. Still structured programming — just better organized.

```fortran
module mod_encapsulated
  implicit none
  private
  public :: gp_model

  type :: gp_model
    real :: lengthscale = 1.0
    real :: noise       = 1e-6
  contains
    procedure :: fit
    procedure :: predict
    procedure :: log_marginal_likelihood
  end type

contains

  subroutine fit(self, X, y)
    class(gp_model), intent(inout) :: self       ! class, not type (required syntax)
    real, intent(in) :: X(:,:), y(:)
  end subroutine

  function predict(self, X_test) result(mu)
    class(gp_model), intent(in) :: self
    real, intent(in) :: X_test(:,:)
    real, allocatable :: mu(:)
  end function

  real function log_marginal_likelihood(self)
    class(gp_model), intent(in) :: self
  end function
end module
```

```fortran
! Usage: clean method-call syntax
type(gp_model) :: model
model%lengthscale = 2.5
call model%fit(X_train, y_train)
mu = model%predict(X_test)
print *, model%log_marginal_likelihood()
```

**Use when:** A type owns behavior that logically belongs to it, but there is only **one variant**. The `class(self)` here is just syntax — no polymorphism happening yet.

---

### Level 2 — OOP: Abstract Type + Extends (Interchangeable Variants)

The abstract type defines the contract. Concrete extensions provide different implementations. Calling code works with `class(abstract_type)` polymorphically.

```fortran
module mod_kernel_oop
  implicit none
  private
  public :: abstract_kernel, rbf_kernel, matern_kernel

  ! ─── Abstract contract ───
  type, abstract :: abstract_kernel
  contains
    procedure(eval_iface), deferred :: eval
    procedure(grad_iface), deferred :: gradient
  end type

  abstract interface
    elemental real function eval_iface(self, x, y)
      import abstract_kernel
      class(abstract_kernel), intent(in) :: self
      real, intent(in) :: x, y
    end function
    elemental real function grad_iface(self, x, y)
      import abstract_kernel
      class(abstract_kernel), intent(in) :: self
      real, intent(in) :: x, y
    end function
  end interface

  ! ─── Concrete variant 1: RBF ───
  type, extends(abstract_kernel) :: rbf_kernel
    real :: lengthscale = 1.0
  contains
    procedure :: eval     => rbf_eval
    procedure :: gradient => rbf_grad
  end type

  ! ─── Concrete variant 2: Matérn 3/2 ───
  type, extends(abstract_kernel) :: matern_kernel
    real :: lengthscale = 1.0
    real :: nu          = 1.5
  contains
    procedure :: eval     => matern_eval
    procedure :: gradient => matern_grad
  end type

contains

  elemental real function rbf_eval(self, x, y)
    class(rbf_kernel), intent(in) :: self
    real, intent(in) :: x, y
    rbf_eval = exp(-0.5 * ((x - y) / self%lengthscale)**2)
  end function

  elemental real function rbf_grad(self, x, y)
    class(rbf_kernel), intent(in) :: self
    real, intent(in) :: x, y
    rbf_grad = -(x - y) / self%lengthscale**2 * exp(-0.5 * ((x - y) / self%lengthscale)**2)
  end function

  elemental real function matern_eval(self, x, y)
    class(matern_kernel), intent(in) :: self
    real, intent(in) :: x, y
    real :: r
    r = sqrt(3.0) * abs(x - y) / self%lengthscale
    matern_eval = (1.0 + r) * exp(-r)
  end function

  elemental real function matern_grad(self, x, y)
    class(matern_kernel), intent(in) :: self
    real, intent(in) :: x, y
    real :: r
    r = sqrt(3.0) * abs(x - y) / self%lengthscale
    matern_grad = -3.0 * (x - y) / self%lengthscale**2 * exp(-r)
  end function
end module
```

```fortran
! ─── The payoff: kernel-agnostic calling code ───
subroutine build_covariance_matrix(kernel, X, K)
  class(abstract_kernel), intent(in) :: kernel   ! polymorphic
  real, intent(in)  :: X(:,:)
  real, intent(out) :: K(:,:)
  integer :: i, j

  do j = 1, size(X, 1)
    do i = 1, size(X, 1)
      K(i,j) = kernel%eval(X(i,1), X(j,1))     ! works for ANY kernel
    end do
  end do
end subroutine
```

```fortran
! Usage: swap kernels without touching build_covariance_matrix
type(rbf_kernel)    :: rbf
type(matern_kernel) :: m32

rbf%lengthscale = 2.5
call build_covariance_matrix(rbf, X, K_rbf)

m32%lengthscale = 3.0
call build_covariance_matrix(m32, X, K_m32)
```

**Use when:** Different types implement the same abstract interface and calling code must work with any of them. Think kernels, optimizers, likelihoods, samplers, preconditioners.

---

### Comparison Table

```
                    ┌─────────────────────────────────────┐
                    │  Does the type represent a CONCEPT   │
                    │  with MULTIPLE interchangeable       │
                    │  implementations?                    │
                    └──────────────┬──────────────────────┘
                                   │
                    ┌──────────────┼──────────────┐
                    │ NO                          │ YES
                    ▼                             ▼
    ┌────────────────────────┐      ┌─────────────────────────┐
    │ Plain derived type or  │      │ Abstract type + extends │
    │ type-bound procedures  │      │ (OOP)                   │
    ├────────────────────────┤      ├─────────────────────────┤
    │ type :: params         │      │ type, abstract :: kernel│
    │   real :: a, b         │      │   procedure(deferred)   │
    │ end type               │      │ end type                │
    │                        │      │                         │
    │ type :: gp             │      │ type, extends(kernel)   │
    │   real :: ls, noise    │      │   :: rbf_kernel         │
    │ contains               │      │   ...                   │
    │   procedure :: fit     │      │ end type                │
    │ end type               │      │                         │
    └────────────────────────┘      └─────────────────────────┘
```

| Question | Plain Type | Type-Bound (no inheritance) | OOP (abstract + extends) |
|---|---|---|---|
| How many variants? | 1 | 1 | 2+ |
| Dispatch | Static (compile-time) | Static (compile-time) | Dynamic (runtime) |
| Self argument | `type(...)` | `class(...)` (syntax only) | `class(...)` (real polymorphism) |
| Example | `type :: point; real :: x,y,z; end type` | `type :: gp; contains; procedure :: fit; end type` | `type, abstract :: optimizer; procedure(deferred) :: step; end type` |
| Overhead | Zero | Zero | Small (vtable lookup) |

---

## 6. Concrete Project Structure Example

```
radio-specmap/
├── fpm.toml
├── src/
│   ├── mod_precision.f90         ! wp = real64, use stdlib kinds
│   ├── mod_kernels.f90           ! elemental kernel functions
│   ├── mod_covariance.f90        ! array: build K matrices from kernels
│   ├── mod_gp.f90                ! structured: GP training & prediction
│   ├── mod_optimizer.f90         ! OOP: abstract_optimizer + lbfgs/adam/...
│   ├── mod_mcmc.f90              ! structured: Metropolis-Hastings, HMC
│   ├── mod_spectrum_map.f90      ! array: spatial interpolation, field reconstruction
│   └── mod_io.f90                ! structured: read/write data
├── test/
│   └── ...
└── app/
    └── main.f90
```

| Module | Paradigm | Rationale |
|---|---|---|
| `mod_kernels` | `elemental` functions | Pure, array-compatible kernel evaluations |
| `mod_covariance`, `mod_gp`, `mod_spectrum_map` | Array + Structured | Matrix math with clear procedural flow |
| `mod_optimizer` | OOP | One abstract type with 3–4 concrete implementations |
| `mod_mcmc` | Structured | Algorithmic pipeline; OOP only if supporting 4+ sampler types |

---

## 7. Summary

| Aspect | Recommendation |
|---|---|
| **Primary paradigm** | **Array + Modular** — write matrix/vector math in array syntax; organize by mathematical concept in modules |
| **Selective OOP** | Only for interchangeable components with 3+ variants (optimizers, kernels, samplers, likelihoods) |
| **Key Fortran features** | `pure`, `elemental`, array intrinsincs, `associate` blocks, `select type` |
| **Paradigms to skip** | Pure functional, deep inheritance hierarchies (>2 levels), heavy generics (until compiler support matures) |
| **How papers map to code** | Papers describe algorithms in structured/array notation → Fortran array syntax is the closest 1:1 mapping |

> **Guiding principle:** let the math drive the paradigm choice. If the paper writes `K = Φ(x)ΣΦ(x)ᵀ`, that's array programming. If it says "we compare three optimizers...", that's OOP. If it describes a step-by-step sampler, that's structured. Don't pick one paradigm and force everything into it — mix them where they fit.
