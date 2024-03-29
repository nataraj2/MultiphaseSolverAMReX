#ifndef AMREX_INTERP_3D_C_H_
#define AMREX_INTERP_3D_C_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_BCRec.H>
#include <AMReX_Vector.H>
#include <cmath>

namespace amrex {

AMREX_GPU_HOST
inline
Vector<Real>
ccinterp_compute_voff (Box const& cbx, IntVect const& ratio, Geometry const& cgeom,
                       Geometry const& fgeom) noexcept
{
    const Box& fbx = amrex::refine(cbx,ratio);
    const auto& flen = amrex::length(fbx);
    const auto& flo = amrex::lbound(fbx);
    const auto& fhi = amrex::ubound(fbx);
    const int nxyz = flen.x + flen.y + flen.z;
    Vector<Real> voff(nxyz);
    
    const Real xrinv = 1./ratio[0];
    const Real yrinv = 1./ratio[1];
    const Real zrinv = 1./ratio[2];

    const Real xtmp = 0.5*(xrinv-1.0);
    const Real ytmp = 0.5*(yrinv-1.0);
    const Real ztmp = 0.5*(zrinv-1.0);

    Real* AMREX_RESTRICT xoff = &voff[0];
    AMREX_PRAGMA_SIMD
    for (int i = flo.x; i <= fhi.x; ++i) {
        const int ic = amrex::coarsen(i, ratio[0]);
        const int ioff = i - ic*ratio[0];
        xoff[i-flo.x] = ioff*xrinv + xtmp;
    }

    Real* AMREX_RESTRICT yoff = &voff[flen.x];
    AMREX_PRAGMA_SIMD
    for (int j = flo.y; j <= fhi.y; ++j) {
        const int jc = amrex::coarsen(j, ratio[1]);
        const int joff = j - jc*ratio[1];
        yoff[j-flo.y] = joff*yrinv + ytmp;
    }

    Real* AMREX_RESTRICT zoff = &voff[flen.x+flen.y];
    AMREX_PRAGMA_SIMD
    for (int k = flo.z; k <= fhi.z; ++k) {
        const int kc = amrex::coarsen(k, ratio[2]);
        const int koff = k - kc*ratio[2];
        zoff[k-flo.z] = koff*zrinv + ztmp;
    }

    return voff;
}

namespace {

AMREX_GPU_HOST_DEVICE AMREX_INLINE void
compute_slopes (const Dim3& lo, const Dim3& hi,
                Array4<Real> const& slopes, const Dim3& slo, const Dim3& shi, int ns,
                Array4<Real const> const& u, int nu,
                BCRec const& bc, int ncomp) noexcept
{
    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                slopes(i,j,k,ns        ) = 0.5*(u(i+1,j,k,nu)-u(i-1,j,k,nu));
                slopes(i,j,k,ns+ncomp  ) = 0.5*(u(i,j+1,k,nu)-u(i,j-1,k,nu));
                slopes(i,j,k,ns+ncomp*2) = 0.5*(u(i,j,k+1,nu)-u(i,j,k-1,nu));
            }
        }
    }

    if (lo.x == slo.x && (bc.lo(0) == BCType::ext_dir || bc.lo(0) == BCType::hoextrap))
    {
        const int i = slo.x;
        if (shi.x-slo.x >= 1) {
            for         (int k = lo.z; k <= hi.z; ++k) {
                for     (int j = lo.y; j <= hi.y; ++j) {
                    slopes(i,j,k,ns) = -(16./15.)*u(i-1,j,k,nu) + 0.5*u(i,j,k,nu)
                        + (2./3.)*u(i+1,j,k,nu) - 0.1*u(i+2,j,k,nu);
                }
            }
        } else {
            for         (int k = lo.z; k <= hi.z; ++k) {
                for     (int j = lo.y; j <= hi.y; ++j) {
                    slopes(i,j,k,ns) = 0.25*(u(i+1,j,k,nu)+5.*u(i,j,k,nu)-6.*u(i-1,j,k,nu));
                }
            }
        }
    }

    if (hi.x == shi.x && (bc.hi(0) == BCType::ext_dir || bc.hi(0) == BCType::hoextrap))
    {
        const int i = shi.x;
        if (shi.x-slo.x >= 1) {
            for         (int k = lo.z; k <= hi.z; ++k) {
                for     (int j = lo.y; j <= hi.y; ++j) {
                    slopes(i,j,k,ns) = (16./15.)*u(i+1,j,k,nu) - 0.5*u(i,j,k,nu)
                        - (2./3.)*u(i-1,j,k,nu) + 0.1*u(i-2,j,k,nu);
                }
            }                
        } else {
            for         (int k = lo.z; k <= hi.z; ++k) {
                for     (int j = lo.y; j <= hi.y; ++j) {
                    slopes(i,j,k,ns) = -0.25*(u(i-1,j,k,nu)+5.*u(i,j,k,nu)-6.*u(i+1,j,k,nu));
                }
            }
        }
    }

    if (lo.y == slo.y && (bc.lo(1) == BCType::ext_dir || bc.lo(1) == BCType::hoextrap))
    {
        const int j = slo.y;
        if (shi.y-slo.y >= 1) {
            for     (int k = lo.z; k <= hi.z; ++k) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    slopes(i,j,k,ns+ncomp) = -(16./15.)*u(i,j-1,k,nu) + 0.5*u(i,j,k,nu)
                        + (2./3.)*u(i,j+1,k,nu) - 0.1*u(i,j+2,k,nu);
                }
            }
        } else {
            for     (int k = lo.z; k <= hi.z; ++k) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    slopes(i,j,k,ns+ncomp) = 0.25*(u(i,j+1,k,nu)+5.*u(i,j,k,nu)-6.*u(i,j-1,k,nu));
                }
            }
        }
    }

    if (hi.y == shi.y && (bc.hi(1) == BCType::ext_dir || bc.hi(1) == BCType::hoextrap))
    {
        const int j = shi.y;
        if (shi.y-slo.y >= 1) {
            for     (int k = lo.z; k <= hi.z; ++k) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    slopes(i,j,k,ns+ncomp) = (16./15.)*u(i,j+1,k,nu) - 0.5*u(i,j,k,nu)
                        - (2./3.)*u(i,j-1,k,nu) + 0.1*u(i,j-2,k,nu);
                }
            }
        } else {
            for     (int k = lo.z; k <= hi.z; ++k) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    slopes(i,j,k,ns+ncomp) = -0.25*(u(i,j-1,k,nu)+5.*u(i,j,k,nu)-6.*u(i,j+1,k,nu));
                }
            }
        }
    }

    if (lo.z == slo.z && (bc.lo(2) == BCType::ext_dir || bc.lo(2) == BCType::hoextrap))
    {
        const int k = slo.z;
        if (shi.z-slo.z >= 1) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    slopes(i,j,k,ns+2*ncomp) = -(16./15.)*u(i,j,k-1,nu) + 0.5*u(i,j,k,nu)
                        + (2./3.)*u(i,j,k+1,nu) - 0.1*u(i,j,k+2,nu);
                }
            }
        } else {
            for     (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    slopes(i,j,k,ns+2*ncomp) = 0.25*(u(i,j,k+1,nu)+5.*u(i,j,k,nu)-6.*u(i,j,k-1,nu));
                }
            }
        }
    }

    if (hi.z == shi.z && (bc.hi(2) == BCType::ext_dir || bc.hi(2) == BCType::hoextrap))
    {
        const int k = shi.z;
        if (shi.z-slo.z >= 1) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    slopes(i,j,k,ns+2*ncomp) = (16./15.)*u(i,j,k+1,nu) - 0.5*u(i,j,k,nu)
                        - (2./3.)*u(i,j,k-1,nu) + 0.1*u(i,j,k-2,nu);
                }
            }
        } else {
            for     (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    slopes(i,j,k,ns+2*ncomp) = -0.25*(u(i,j,k-1,nu)+5.*u(i,j,k,nu)-6.*u(i,j,k+1,nu));
                }
            }
        }
    }
}

}

AMREX_GPU_HOST_DEVICE inline void
cellconslin_slopes_linlim (Box const& bx, Array4<Real> const& slopes,
                           Array4<Real const> const& u, const int icomp, const int ncomp,
                           BCRec const* AMREX_RESTRICT bcr) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    const auto slo  = amrex::lbound(slopes);
    const auto shi  = amrex::ubound(slopes);

    const Array4<Real> sf(slopes, ncomp*AMREX_SPACEDIM);  // slope factor

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                sf(i,j,k,0) = 1.0;
                sf(i,j,k,1) = 1.0;
                sf(i,j,k,2) = 1.0;
            }
        }
    }

    for (int n = 0; n < ncomp; ++n)
    {
        const int nu = n + icomp;
        compute_slopes(lo, hi, slopes, slo, shi, n, u, nu, bcr[n], ncomp);

        for         (int k = lo.z; k <= hi.z; ++k) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    Real cen  = slopes(i,j,k,n);
                    Real forw = 2.0*(u(i+1,j,k,nu)-u(i  ,j,k,nu));
                    Real back = 2.0*(u(i  ,j,k,nu)-u(i-1,j,k,nu));
                    Real slp = (forw*back >= 0.0) ? amrex::min(std::abs(forw),std::abs(back)) : 0.0;
                    slopes(i,j,k,n) = std::copysign(1.0,cen)*amrex::min(slp,std::abs(cen));
                    if (cen != 0.0) {
                        sf(i,j,k,0) = amrex::min(sf(i,j,k,0), slopes(i,j,k,n)/cen);
                    } else {
                        sf(i,j,k,0) = 0.0;
                    }

                    cen  = slopes(i,j,k,n+ncomp);
                    forw = 2.0*(u(i,j+1,k,nu)-u(i,j  ,k,nu));
                    back = 2.0*(u(i,j  ,k,nu)-u(i,j-1,k,nu));
                    slp = (forw*back >= 0.0) ? amrex::min(std::abs(forw),std::abs(back)) : 0.0;
                    slopes(i,j,k,n+ncomp) = std::copysign(1.0,cen)*amrex::min(slp,std::abs(cen));
                    if (cen != 0.0) {
                        sf(i,j,k,1) = amrex::min(sf(i,j,k,1), slopes(i,j,k,n+ncomp)/cen);
                    } else {
                        sf(i,j,k,1) = 0.0;
                    }

                    cen  = slopes(i,j,k,n+ncomp*2);
                    forw = 2.0*(u(i,j,k+1,nu)-u(i,j,k  ,nu));
                    back = 2.0*(u(i,j,k  ,nu)-u(i,j,k-1,nu));
                    slp = (forw*back >= 0.0) ? amrex::min(std::abs(forw),std::abs(back)) : 0.0;
                    slopes(i,j,k,n+ncomp*2) = std::copysign(1.0,cen)*amrex::min(slp,std::abs(cen));
                    if (cen != 0.0) {
                        sf(i,j,k,2) = amrex::min(sf(i,j,k,2), slopes(i,j,k,n+ncomp*2)/cen);
                    } else {
                        sf(i,j,k,2) = 0.0;
                    }
                }
            }
        }
    }

    for             (int n = 0; n < ncomp; ++n) {
        for         (int k = lo.z; k <= hi.z; ++k) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    slopes(i,j,k,n        ) *= sf(i,j,k,0);
                    slopes(i,j,k,n+ncomp  ) *= sf(i,j,k,1);
                    slopes(i,j,k,n+ncomp*2) *= sf(i,j,k,2);
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE inline void
cellconslin_interp (Box const& bx,
                    Array4<Real> const& fine, const int fcomp, const int ncomp,
                    Array4<Real const> const& slopes,
                    Array4<Real const> const& crse, const int ccomp,
                    Real const* AMREX_RESTRICT voff, IntVect const& ratio) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    Box vbox(slopes);
    vbox.refine(ratio);
    const auto vlo  = amrex::lbound(vbox);
    const auto vlen = amrex::length(vbox);
    Real const* AMREX_RESTRICT xoff = voff;
    Real const* AMREX_RESTRICT yoff = voff + vlen.x;
    Real const* AMREX_RESTRICT zoff = voff + (vlen.x+vlen.y);

    for (int n = 0; n < ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
            const int kc = amrex::coarsen(k,ratio[2]);
            for (int j = lo.y; j <= hi.y; ++j) {
                const int jc = amrex::coarsen(j,ratio[1]);
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    const int ic = amrex::coarsen(i,ratio[0]);
                    fine(i,j,k,n+fcomp) = crse(ic,jc,kc,n+ccomp)
                        + xoff[i-vlo.x] * slopes(ic,jc,kc,n)
                        + yoff[j-vlo.y] * slopes(ic,jc,kc,n+ncomp)
                        + zoff[k-vlo.z] * slopes(ic,jc,kc,n+2*ncomp);
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE inline void
cellconslin_slopes_mclim (Box const& bx, Array4<Real> const& slopes,
                          Array4<Real const> const& u, const int icomp, const int ncomp,
                          BCRec const* AMREX_RESTRICT bcr) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    const auto slo = amrex::lbound(slopes);
    const auto shi = amrex::ubound(slopes);

    const Array4<Real> mm(slopes, ncomp*AMREX_SPACEDIM);  // min and max

    for             (int n = 0; n < ncomp; ++n) {
        const int nu = n + icomp;
        for         (int k = lo.z; k <= hi.z; ++k) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    Real cmn = u(i,j,k,nu);
                    Real cmx = cmn;
                    for         (int koff = -1; koff <= 1; ++koff) {
                        for     (int joff = -1; joff <= 1; ++joff) {
                            for (int ioff = -1; ioff <= 1; ++ioff) {
                                cmn = amrex::min(cmn,u(i+ioff,j+joff,k+koff,nu));
                                cmx = amrex::max(cmx,u(i+ioff,j+joff,k+koff,nu));
                            }
                        }
                    }
                    mm(i,j,k,n      ) = cmn - u(i,j,k,nu);
                    mm(i,j,k,n+ncomp) = cmx - u(i,j,k,nu);
                }
            }
        }

        compute_slopes(lo, hi, slopes, slo, shi, n, u, nu, bcr[n], ncomp);

        for         (int k = lo.z; k <= hi.z; ++k) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    Real cen  = slopes(i,j,k,n);
                    Real forw = 2.0*(u(i+1,j,k,nu)-u(i  ,j,k,nu));
                    Real back = 2.0*(u(i  ,j,k,nu)-u(i-1,j,k,nu));
                    Real slp = (forw*back >= 0.0) ? amrex::min(std::abs(forw),std::abs(back)) : 0.0;
                    slopes(i,j,k,n) = std::copysign(1.0,cen)*amrex::min(slp,std::abs(cen));

                    cen  = slopes(i,j,k,n+ncomp);
                    forw = 2.0*(u(i,j+1,k,nu)-u(i,j  ,k,nu));
                    back = 2.0*(u(i,j  ,k,nu)-u(i,j-1,k,nu));
                    slp = (forw*back >= 0.0) ? amrex::min(std::abs(forw),std::abs(back)) : 0.0;
                    slopes(i,j,k,n+ncomp) = std::copysign(1.0,cen)*amrex::min(slp,std::abs(cen));

                    cen  = slopes(i,j,k,n+ncomp*2);
                    forw = 2.0*(u(i,j,k+1,nu)-u(i,j,k  ,nu));
                    back = 2.0*(u(i,j,k  ,nu)-u(i,j,k-1,nu));
                    slp = (forw*back >= 0.0) ? amrex::min(std::abs(forw),std::abs(back)) : 0.0;
                    slopes(i,j,k,n+ncomp*2) = std::copysign(1.0,cen)*amrex::min(slp,std::abs(cen));
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE inline void
cellconslin_fine_alpha (Box const& bx, Array4<Real> const& alpha,
                        Array4<Real const> const& slopes, const int ncomp,
                        Real const* AMREX_RESTRICT voff, IntVect const& ratio) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    const auto slo = amrex::lbound(slopes);
    const auto shi = amrex::ubound(slopes);

    const Array4<Real const> mm(slopes, ncomp*AMREX_SPACEDIM);  // min and max

    const auto vlo  = amrex::lbound(alpha);
    const auto vlen = amrex::length(alpha);
    Real const* AMREX_RESTRICT xoff = voff;
    Real const* AMREX_RESTRICT yoff = voff + vlen.x;
    Real const* AMREX_RESTRICT zoff = voff + (vlen.x+vlen.y);

    for (int n = 0; n < ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
            const int kc = amrex::coarsen(k,ratio[2]);
            for (int j = lo.y; j <= hi.y; ++j) {
                const int jc = amrex::coarsen(j,ratio[1]);
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    const int ic = amrex::coarsen(i,ratio[0]);
                    const Real dummy_fine = xoff[i-vlo.x]*slopes(ic,jc,kc,n)
                        +                   yoff[j-vlo.y]*slopes(ic,jc,kc,n+ncomp)
                        +                   zoff[k-vlo.z]*slopes(ic,jc,kc,n+ncomp*2);

                    if (dummy_fine > mm(ic,jc,kc,n+ncomp) && dummy_fine != 0.0) {
                        alpha(i,j,k,n) = mm(ic,jc,kc,n+ncomp) / dummy_fine;
                    } else if (dummy_fine < mm(ic,jc,kc,n) && dummy_fine != 0.0) {
                        alpha(i,j,k,n) = mm(ic,jc,kc,n) / dummy_fine;
                    } else {
                        alpha(i,j,k,n) = 1.0;
                    }
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE inline void
cellconslin_slopes_mmlim (Box const& bx, Array4<Real> const& slopes,
                          Array4<Real const> const& alpha,
                          const int ncomp, IntVect const& ratio) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    for (int n = 0; n < ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
            const int kk = k*ratio[2];
            for (int j = lo.y; j <= hi.y; ++j) {
                const int jj = j*ratio[1];
                for (int i = lo.x; i <= hi.x; ++i) {
                    const int ii = i*ratio[0];
                    Real a = 1.0;
                    for         (int koff = 0; koff < ratio[2]; ++koff) {
                        for     (int joff = 0; joff < ratio[1]; ++joff) {
                            for (int ioff = 0; ioff < ratio[0]; ++ioff) {
                                a = amrex::min(a, alpha(ii+ioff,jj+joff,kk+koff,n));
                            }
                        }
                    }
                    slopes(i,j,k,n        ) *= a;
                    slopes(i,j,k,n+ncomp  ) *= a;
                    slopes(i,j,k,n+ncomp*2) *= a;
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE inline void
pcinterp_interp (Box const& bx,
                 Array4<Real> const& fine, const int fcomp, const int ncomp,
                 Array4<Real const> const& crse, const int ccomp,
                 IntVect const& ratio) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    for (int n = 0; n < ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
            const int kc = amrex::coarsen(k,ratio[2]);
            for (int j = lo.y; j <= hi.y; ++j) {
                const int jc = amrex::coarsen(j,ratio[1]);
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    const int ic = amrex::coarsen(i,ratio[0]);
                    fine(i,j,k,n+fcomp) = crse(ic,jc,kc,n+ccomp);
                }
            }
        }
    }
}

namespace {
    static constexpr int ix   = 0;
    static constexpr int iy   = 1;
    static constexpr int iz   = 2;
    static constexpr int ixy  = 3;
    static constexpr int ixz  = 4;
    static constexpr int iyz  = 5;
    static constexpr int ixyz = 6;
}

AMREX_GPU_HOST_DEVICE inline void
nodebilin_slopes (Box const& bx, Array4<Real> const& slope, Array4<Real const> const& u,
                  const int icomp, const int ncomp, IntVect const& ratio) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    const Real rx = 1.0/ratio[0];
    const Real ry = 1.0/ratio[1];
    const Real rz = 1.0/ratio[2];

    for (int n = 0; n < ncomp; ++n) {
        const int nu = n + icomp;
        for         (int k = lo.z; k <= hi.z; ++k) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    Real dx00 = u(i+1,j,k,nu) - u(i,j,k,nu);
                    Real d0x0 = u(i,j+1,k,nu) - u(i,j,k,nu);
                    Real d00x = u(i,j,k+1,nu) - u(i,j,k,nu);

                    Real dx10 = u(i+1,j+1,k,nu) - u(i,j+1,k,nu);
                    Real dx01 = u(i+1,j,k+1,nu) - u(i,j,k+1,nu);
                    Real d0x1 = u(i,j+1,k+1,nu) - u(i,j,k+1,nu);

                    Real dx11 = u(i+1,j+1,k+1,nu) - u(i,j+1,k+1,nu);

                    slope(i,j,k,n+ncomp*ix  ) = rx*dx00;
                    slope(i,j,k,n+ncomp*iy  ) = ry*d0x0;
                    slope(i,j,k,n+ncomp*iz  ) = rz*d00x;
                    slope(i,j,k,n+ncomp*ixy ) = rx*ry*(dx10 - dx00);
                    slope(i,j,k,n+ncomp*ixz ) = rx*rz*(dx01 - dx00);
                    slope(i,j,k,n+ncomp*iyz ) = ry*rz*(d0x1 - d0x0);
                    slope(i,j,k,n+ncomp*ixyz) = rx*ry*rz*(dx11 - dx01 - dx10 + dx00);
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE inline void
nodebilin_interp (Box const& bx, Array4<Real> const& fine, const int fcomp, const int ncomp,
                  Array4<Real const> const& slope, Array4<Real const> const& crse,
                  const int ccomp, IntVect const& ratio) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    const auto chi = amrex::ubound(slope);

    for (int n = 0; n < ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
            const int kc = amrex::min(amrex::coarsen(k,ratio[2]),chi.z);
            const Real fz = k - kc*ratio[2];
            for (int j = lo.y; j <= hi.y; ++j) {
                const int jc = amrex::min(amrex::coarsen(j,ratio[1]),chi.y);
                const Real fy = j - jc*ratio[1];
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    const int ic = amrex::min(amrex::coarsen(i,ratio[0]),chi.x);
                    const Real fx = i - ic*ratio[0];
                    fine(i,j,k,n+fcomp) = crse(ic,jc,kc,n+ccomp)
                        + fx*slope(ic,jc,kc,n+ncomp*ix)
                        + fy*slope(ic,jc,kc,n+ncomp*iy)
                        + fz*slope(ic,jc,kc,n+ncomp*iz)
                        + fx*fy*slope(ic,jc,kc,n+ncomp*ixy)
                        + fx*fz*slope(ic,jc,kc,n+ncomp*ixz)
                        + fy*fz*slope(ic,jc,kc,n+ncomp*iyz)
                        + fx*fy*fz*slope(ic,jc,kc,n+ncomp*ixyz);
                }
            }
        }
    }
}

}

#endif
