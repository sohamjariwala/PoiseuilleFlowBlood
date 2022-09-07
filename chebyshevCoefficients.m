function x = chebyshevCoefficients ( a, b, n )
angle = ( 1 : 2 : ( 2 * n - 1 ) ) * pi / ( 2 * n );
angle = angle';
x = -cos ( angle );
x = 0.5 * ( a + b ) + x * 0.5 * ( b - a );
end