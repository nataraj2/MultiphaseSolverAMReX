! ========================================= !
! Look-up tables for computational geometry !
! ========================================= !
module compgeom_lookup
  implicit none
  
  ! Number of new vertices on cut plane 
  integer, dimension(8) :: cutTri_nvert
  data cutTri_nvert(1:8) / 0, 2, 2, 2, 2, 2, 2, 0/
  ! Number of resulting tris on positive side 
  integer, dimension(8) :: cutTri_np  
  data cutTri_np(1:8)   / 0, 1, 1, 2, 1, 2, 2, 1/
   ! Number of resulting tris on negative side
  integer, dimension(8) :: cutTri_nn  
  data cutTri_nn(1:8)   / 1, 2, 2, 1, 2, 1, 1, 0/
  ! First point on intersection 
  integer, dimension(2,8) :: cutTri_v1
  data cutTri_v1(1:2,1) /-1,-1/
  data cutTri_v1(1:2,2) / 1, 1/
  data cutTri_v1(1:2,3) / 2, 2/
  data cutTri_v1(1:2,4) / 3, 3/
  data cutTri_v1(1:2,5) / 3, 3/ 
  data cutTri_v1(1:2,6) / 2, 2/
  data cutTri_v1(1:2,7) / 1, 1/
  data cutTri_v1(1:2,8) /-1,-1/
  ! Second point on intersection 
  integer, dimension(2,8) :: cutTri_v2
  data cutTri_v2(1:2,1) /-1,-1/
  data cutTri_v2(1:2,2) / 2, 3/
  data cutTri_v2(1:2,3) / 1, 3/
  data cutTri_v2(1:2,4) / 1, 2/
  data cutTri_v2(1:2,5) / 1, 2/ 
  data cutTri_v2(1:2,6) / 1, 3/
  data cutTri_v2(1:2,7) / 2, 3/
  data cutTri_v2(1:2,8) /-1,-1/
  ! Vertices in each tri
  integer, dimension(3,3,8) :: cutTri_v
  data cutTri_v(1:3,1:3,1) / 1, 2, 3, -1,-1,-1, -1,-1,-1/
  data cutTri_v(1:3,1:3,2) / 1, 4, 5,  2, 5, 4,  2, 3, 5/
  data cutTri_v(1:3,1:3,3) / 2, 5, 4,  1, 4, 5,  1, 5, 3/
  data cutTri_v(1:3,1:3,4) / 1, 2, 5,  1, 5, 4,  3, 4, 5/
  data cutTri_v(1:3,1:3,5) / 3, 4, 5,  1, 2, 5,  1, 5, 4/
  data cutTri_v(1:3,1:3,6) / 1, 4, 5,  1, 5, 3,  2, 5, 4/
  data cutTri_v(1:3,1:3,7) / 2, 3, 5,  2, 5, 4,  1, 4, 5/
  data cutTri_v(1:3,1:3,8) / 1, 2, 3, -1,-1,-1, -1,-1,-1/
  
  ! Cell vertices to edges
  integer, dimension(2,12) :: verts2edge
  data verts2edge(:, 1) / 1, 2 /
  data verts2edge(:, 2) / 1, 3 /
  data verts2edge(:, 3) / 1, 5 /
  data verts2edge(:, 4) / 2, 4 /
  data verts2edge(:, 5) / 2, 6 /
  data verts2edge(:, 6) / 3, 4 /
  data verts2edge(:, 7) / 3, 7 /
  data verts2edge(:, 8) / 4, 8 /
  data verts2edge(:, 9) / 5, 6 /
  data verts2edge(:,10) / 5, 7 /
  data verts2edge(:,11) / 6, 8 /
  data verts2edge(:,12) / 7, 8 /
  
  ! Vertices on cell to 5 tets
  integer, dimension(4,5) :: verts2tets
  data verts2tets(:,1) / 1, 4, 3, 7 /
  data verts2tets(:,2) / 1, 2, 4, 6 /
  data verts2tets(:,3) / 1, 5, 6, 7 /
  data verts2tets(:,4) / 1, 6, 4, 7 /
  data verts2tets(:,5) / 4, 7, 6, 8 /
  
  ! Vertices on cell to 8 tets for flux-based SL
  integer, dimension(4,8) :: tet_map
  data tet_map(1:4,1) / 7, 4, 3, 6 /
  data tet_map(1:4,2) / 6, 3, 2, 4 /
  data tet_map(1:4,3) / 6, 2, 1, 4 /
  data tet_map(1:4,4) / 7, 8, 4, 6 /
  data tet_map(1:4,5) / 6, 5, 8, 4 /
  data tet_map(1:4,6) / 6, 5, 4, 1 /
  data tet_map(1:4,7) / 5, 6, 8, 9 /
  data tet_map(1:4,8) / 6, 7, 8, 9 /

  ! Vertices on cell to 6 tets for pure semi-Lagrangian
  integer, dimension(4,6) :: tet_map2
  data tet_map2(1:4,1) / 7, 4, 3, 6 /
  data tet_map2(1:4,2) / 6, 3, 2, 4 /
  data tet_map2(1:4,3) / 6, 2, 1, 4 /
  data tet_map2(1:4,4) / 7, 8, 4, 6 /
  data tet_map2(1:4,5) / 6, 5, 8, 4 /
  data tet_map2(1:4,6) / 6, 5, 4, 1 /

  ! Vertices on cell to 8 tets for flux-based SL
  ! x(i) face
  integer, dimension(4,8), parameter :: tet_mapxm=reshape((/ &
  7, 16, 15,  6,&
  6, 15, 14, 16,&
  6, 14, 13, 16,&
  7,  8, 16,  6,&
  6,  5,  8, 16,&
  6,  5, 16, 13,&
  5,  6,  8, 17,&
  6,  7,  8, 17/),(/4,8/))

  ! Vertices on cell to 8 tets for flux-based SL
  ! x(i+1) face
  integer, dimension(4,8), parameter :: tet_mapxp=reshape((/ &
  3, 12, 11,  2,&
  2, 11, 10, 12,&
  2, 10,  9, 12,&
  3,  4, 12,  2,&
  2,  1,  4, 12,&
  2,  1, 12,  9,&
  2,  3,  4, 17,&
  1,  2,  4, 17/),(/4,8/))

  ! Vertices on cell to 8 tets for flux-based SL
  ! y(j) face
  integer, dimension(4,8), parameter :: tet_mapym=reshape((/ &
  2,  9, 10,  6,&
  6, 10, 14,  9,&
  6, 14, 13,  9,&
  2,  1,  9,  6,&
  6,  5,  1,  9,&
  6,  5,  9, 13,&
  5,  6,  1, 17,&
  6,  2,  1, 17/),(/4,8/))

  ! Vertices on cell to 8 tets for flux-based SL
  ! y(j+1) face
  integer, dimension(4,8), parameter :: tet_mapyp=reshape((/ &
  3, 12, 11,  7,&
  7, 11, 15, 12,&
  7, 15, 16, 12,&
  3,  4, 12,  7,&
  7,  8,  4, 12,&
  7,  8, 12, 16,&
  8,  7,  4, 17,&
  7,  3,  4, 17/),(/4,8/))

  ! Vertices on cell to 8 tets for flux-based SL
  ! z(k) face
  integer, dimension(4,8), parameter :: tet_mapzm=reshape((/ &
  8, 12, 16,  5,&
  5, 16, 13, 12,&
  5, 13,  9, 12,&
  8,  4, 12,  5,&
  5,  1,  4, 12,&
  5,  1, 12,  9,&
  1,  5,  4, 17,&
  5,  8,  4, 17/),(/4,8/))

  ! Vertices on cell to 8 tets for flux-based SL
  ! z(k+1) face
  integer, dimension(4,8), parameter :: tet_mapzp=reshape((/ &
  7, 11, 15,  6,&
  6, 15, 14, 11,&
  6, 14, 10, 11,&
  7,  3, 11,  6,&
  6,  2,  3, 11,&
  6,  2, 11, 10,&
  2,  6,  3, 17,&
  6,  7,  3, 17/),(/4,8/))

!##################
! Vertices on cell to 8 tets for flux-based SL
! x(i) face
integer, dimension(4,2), parameter :: tet_corxm=reshape((/ &
5,  6,  8, 17,&
6,  7,  8, 17/),(/4,2/))

! Vertices on cell to 8 tets for flux-based SL
! x(i+1) face
integer, dimension(4,2), parameter :: tet_corxp=reshape((/ &
2,  3,  4, 17,&
1,  2,  4, 17/),(/4,2/))

! Vertices on cell to 8 tets for flux-based SL
! y(j) face
integer, dimension(4,2), parameter :: tet_corym=reshape((/ &
5,  6,  1, 17,&
6,  2,  1, 17/),(/4,2/))

! Vertices on cell to 8 tets for flux-based SL
! y(j+1) face
integer, dimension(4,2), parameter :: tet_coryp=reshape((/ &
8,  7,  4, 17,&
7,  3,  4, 17/),(/4,2/))

! Vertices on cell to 8 tets for flux-based SL
! z(k) face
integer, dimension(4,2), parameter :: tet_corzm=reshape((/ &
1,  5,  4, 17,&
5,  8,  4, 17/),(/4,2/))

! Vertices on cell to 8 tets for flux-based SL
! z(k+1) face
integer, dimension(4,2), parameter :: tet_corzp=reshape((/ &
2,  6,  3, 17,&
6,  7,  3, 17/),(/4,2/))
!####################

  ! Vertices on cell face to 2 tris
  integer, dimension(3,2) :: verts2tris
  data verts2tris(:,1) / 1, 2, 4 /
  data verts2tris(:,2) / 1, 4, 3 /
  
  ! Number of new vertices on cut plane
  integer, dimension(16) :: cut_nvert
  data cut_nvert(1:16) / 0, 3, 3, 4, 3, 4, 4, 3, 3, 4, 4, 3, 4, 3, 3, 0/
  ! Number of resulting tets
  integer, dimension(16) :: cut_ntets
  data cut_ntets(1:16) / 1, 4, 4, 6, 4, 6, 6, 4, 4, 6, 6, 4, 6, 4, 4, 1/

  ! First point on intersection 
  integer, dimension(4,16) :: cut_v1
  data cut_v1(1:4, 1) /-1,-1,-1,-1/
  data cut_v1(1:4, 2) / 1, 1, 1,-1/
  data cut_v1(1:4, 3) / 2, 2, 2,-1/
  data cut_v1(1:4, 4) / 1, 2, 1, 2/
  data cut_v1(1:4, 5) / 3, 3, 3,-1/
  data cut_v1(1:4, 6) / 1, 3, 1, 3/
  data cut_v1(1:4, 7) / 2, 3, 2, 3/
  data cut_v1(1:4, 8) / 4, 4, 4,-1/
  data cut_v1(1:4, 9) / 4, 4, 4,-1/
  data cut_v1(1:4,10) / 1, 4, 1, 4/
  data cut_v1(1:4,11) / 2, 4, 2, 4/
  data cut_v1(1:4,12) / 3, 3, 3,-1/
  data cut_v1(1:4,13) / 3, 4, 3, 4/
  data cut_v1(1:4,14) / 2, 2, 2,-1/
  data cut_v1(1:4,15) / 1, 1, 1,-1/
  data cut_v1(1:4,16) /-1,-1,-1,-1/

  ! Second point on intersection 
  integer, dimension(4,16) :: cut_v2
  data cut_v2(1:4, 1) /-1,-1,-1,-1/
  data cut_v2(1:4, 2) / 2, 3, 4,-1/
  data cut_v2(1:4, 3) / 3, 4, 1,-1/
  data cut_v2(1:4, 4) / 4, 4, 3, 3/
  data cut_v2(1:4, 5) / 4, 1, 2,-1/
  data cut_v2(1:4, 6) / 4, 4, 2, 2/
  data cut_v2(1:4, 7) / 4, 4, 1, 1/
  data cut_v2(1:4, 8) / 1, 2, 3,-1/
  data cut_v2(1:4, 9) / 1, 2, 3,-1/
  data cut_v2(1:4,10) / 3, 3, 2, 2/
  data cut_v2(1:4,11) / 3, 3, 1, 1/
  data cut_v2(1:4,12) / 4, 1, 2,-1/
  data cut_v2(1:4,13) / 2, 2, 1, 1/
  data cut_v2(1:4,14) / 3, 4, 1,-1/
  data cut_v2(1:4,15) / 2, 3, 4,-1/
  data cut_v2(1:4,16) /-1,-1,-1,-1/
  
  ! Side of cut plane (used to update i,j,k)
  integer, dimension(6,  16) :: cut_side
  data cut_side(1:6, 1) / 1,-1,-1,-1,-1,-1/
  data cut_side(1:6, 2) / 2, 1, 1, 1,-1,-1/
  data cut_side(1:6, 3) / 2, 1, 1, 1,-1,-1/
  data cut_side(1:6, 4) / 2, 2, 2, 1, 1, 1/
  data cut_side(1:6, 5) / 2, 1, 1, 1,-1,-1/
  data cut_side(1:6, 6) / 2, 2, 2, 1, 1, 1/
  data cut_side(1:6, 7) / 2, 2, 2, 1, 1, 1/
  data cut_side(1:6, 8) / 2, 2, 2, 1,-1,-1/
  data cut_side(1:6, 9) / 2, 1, 1, 1,-1,-1/
  data cut_side(1:6,10) / 2, 2, 2, 1, 1, 1/
  data cut_side(1:6,11) / 2, 2, 2, 1, 1, 1/
  data cut_side(1:6,12) / 2, 2, 2, 1,-1,-1/
  data cut_side(1:6,13) / 2, 2, 2, 1, 1, 1/
  data cut_side(1:6,14) / 2, 2, 2, 1,-1,-1/
  data cut_side(1:6,15) / 2, 2, 2, 1,-1,-1/
  data cut_side(1:6,16) / 2,-1,-1,-1,-1,-1/

  ! Vertices in each tet
  integer, dimension(4,6,16) :: cut_vtet
  data cut_vtet(1:4,1:6, 1) / 1, 2, 3, 4, -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1/
  data cut_vtet(1:4,1:6, 2) / 5, 7, 6, 1,  6, 2, 3, 4,  4, 2, 5, 6,  5, 6, 7, 4, -1,-1,-1,-1, -1,-1,-1,-1/
  data cut_vtet(1:4,1:6, 3) / 7, 5, 6, 2,  1, 3, 4, 6,  1, 5, 3, 6,  5, 7, 6, 1, -1,-1,-1,-1, -1,-1,-1,-1/
  data cut_vtet(1:4,1:6, 4) / 5, 8, 6, 2,  5, 7, 8, 1,  5, 1, 8, 2,  5, 6, 8, 4,  5, 8, 7, 3,  5, 8, 3, 4/
  data cut_vtet(1:4,1:6, 5) / 6, 5, 7, 3,  2, 1, 4, 6,  6, 5, 4, 2,  6, 7, 5, 2, -1,-1,-1,-1, -1,-1,-1,-1/
  data cut_vtet(1:4,1:6, 6) / 5, 6, 8, 3,  5, 8, 7, 1,  5, 8, 1, 3,  5, 8, 6, 4,  5, 7, 8, 2,  5, 8, 4, 2/
  data cut_vtet(1:4,1:6, 7) / 8, 6, 5, 3,  5, 7, 8, 2,  8, 5, 2, 3,  8, 5, 6, 4,  5, 8, 7, 1,  5, 8, 1, 4/
  data cut_vtet(1:4,1:6, 8) / 1, 2, 3, 7,  1, 2, 7, 6,  5, 7, 6, 1,  5, 6, 7, 4, -1,-1,-1,-1, -1,-1,-1,-1/
  data cut_vtet(1:4,1:6, 9) / 5, 6, 7, 4,  1, 2, 3, 6,  5, 1, 3, 6,  5, 7, 6, 3, -1,-1,-1,-1, -1,-1,-1,-1/
  data cut_vtet(1:4,1:6,10) / 5, 8, 6, 4,  5, 7, 8, 1,  5, 8, 4, 1,  5, 6, 8, 3,  5, 8, 7, 2,  5, 8, 2, 3/
  data cut_vtet(1:4,1:6,11) / 8, 5, 6, 4,  5, 8, 7, 2,  8, 2, 5, 4,  8, 6, 5, 3,  5, 7, 8, 1,  5, 8, 3, 1/
  data cut_vtet(1:4,1:6,12) / 1, 4, 2, 7,  4, 1, 6, 7,  6, 7, 5, 4,  6, 5, 7, 3, -1,-1,-1,-1, -1,-1,-1,-1/
  data cut_vtet(1:4,1:6,13) / 8, 6, 5, 4,  5, 7, 8, 3,  8, 4, 5, 3,  8, 5, 6, 2,  5, 8, 7, 1,  5, 8, 1, 2/
  data cut_vtet(1:4,1:6,14) / 3, 4, 1, 7,  7, 6, 3, 4,  7, 6, 5, 3,  7, 5, 6, 2, -1,-1,-1,-1, -1,-1,-1,-1/
  data cut_vtet(1:4,1:6,15) / 7, 4, 2, 3,  2, 3, 6, 7,  5, 6, 7, 2,  5, 7, 6, 1, -1,-1,-1,-1, -1,-1,-1,-1/
  data cut_vtet(1:4,1:6,16) / 1, 2, 3, 4, -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1/

  ! Number of triangles on cut plane
  integer, dimension(16) :: cut_ntris
  data cut_ntris(1:16) / 0, 1, 1, 2, 1, 2, 2, 1, 1, 2, 2, 1, 2, 1, 1, 0/
  
  ! Vertices in each tri on cut plane
  integer, dimension(3,2,16) :: cut_vtri
  data cut_vtri(1:3,1:2, 1) /-1,-1,-1, -1,-1,-1/
  data cut_vtri(1:3,1:2, 2) / 5, 7, 6, -1,-1,-1/
  data cut_vtri(1:3,1:2, 3) / 5, 6, 7, -1,-1,-1/
  data cut_vtri(1:3,1:2, 4) / 5, 8, 6,  5, 7, 8/
  data cut_vtri(1:3,1:2, 5) / 5, 7, 6, -1,-1,-1/
  data cut_vtri(1:3,1:2, 6) / 5, 6, 8,  5, 8, 7/
  data cut_vtri(1:3,1:2, 7) / 5, 8, 6,  5, 7, 8/
  data cut_vtri(1:3,1:2, 8) / 5, 7, 6, -1,-1,-1/
  data cut_vtri(1:3,1:2, 9) / 5, 6, 7, -1,-1,-1/
  data cut_vtri(1:3,1:2,10) / 5, 8, 6,  5, 7, 8/
  data cut_vtri(1:3,1:2,11) / 5, 6, 8,  5, 8, 7/
  data cut_vtri(1:3,1:2,12) / 5, 6, 7, -1,-1,-1/
  data cut_vtri(1:3,1:2,13) / 5, 8, 6,  5, 7, 8/
  data cut_vtri(1:3,1:2,14) / 5, 7, 6, -1,-1,-1/
  data cut_vtri(1:3,1:2,15) / 5, 6, 7, -1,-1,-1/
  data cut_vtri(1:3,1:2,16) /-1,-1,-1, -1,-1,-1/
  
  ! Index of first positive tet = # tets - # negative tets + 1
  integer, dimension(16) :: cut_nntet
  data cut_nntet(1:16) / 1, 2, 2, 4, 2, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 2/
  
end module compgeom_lookup
