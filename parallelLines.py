import math


def point_between_parallel_lines(x0, y0, x1, y1, x2, y2, xp, yp):
    A = (x1 - x0, y1 - y0)  # Vector A

    # Vector from the first line to the point
    B1 = (xp - x0, yp - y0)
    # Vector from the second line to the point
    B2 = (xp - x2, yp - y2)

    # Z-component of the cross product for the first line
    z1 = A[0] * B1[1] - A[1] * B1[0]
    # Z-component of the cross product for the second line
    z2 = A[0] * B2[1] - A[1] * B2[0]

    # If z1 and z2 have different signs, the point is outside
    if (z1 * z2 < 0):
        return False
    else:
        return True
    return
'''
# Example usage:x0, y0 = 0, 0  
x0 = 9,y0 = 5
print(x0)
x1 = 12
y1 = 6

x2 = 2
y2 = 5

xp = 3
yp = 7
'''
# Define the coordinates of the first line
x0, y0 = 0, 0  # Point 1 of the first line
x1, y1 = 10, 0  # Point 2 of the first line

# Define the coordinates of a point on the second line
# Assuming the second line is parallel and offset in the y direction by some amount
y_offset = 5
x2, y2 = x0, y0 + y_offset  # Point on the second line with the same x but an offset y

# Define the coordinates of the point to check
xp, yp = 5, 6  # The point to check

is_between = point_between_parallel_lines(x0, y0, x1, y1, x2, y2, xp, yp)
print("The point is between the parallel lines." if is_between else "The point is not between the parallel lines.")

