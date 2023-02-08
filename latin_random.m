function x = latin_random ( dim_num, point_num )

%*****************************************************************************80
%
%% latin_random() returns points in a Latin Random Square.
%
%  Discussion:
%
%    In each spatial dimension, there will be exactly one
%    point whose coordinate value lies between consecutive
%    values in the list:
%
%      ( 0, 1, 2, ..., point_num ) / point_num
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    09 July 2021
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    integer DIM_NUM, the spatial dimension.
%
%    integer POINT_NUM, the number of points.
%
%  Output:
%
%    real X(DIM_NUM,POINT_NUM), the points.
%
  x = rand ( dim_num, point_num );
%
%  For spatial dimension I,
%    pick a random permutation of 1 to POINT_NUM,
%    force the corresponding I-th components of X to lie in the
%    interval ( PERM(J)-1, PERM(J) ) / POINT_NUM.
%
  for i = 1: dim_num

    perm = randperm ( point_num );

    for j = 1: point_num
      x(i,j) = ( ( perm(j) - 1 ) + x(i,j) ) / ( point_num );
    end

  end

  return
end
