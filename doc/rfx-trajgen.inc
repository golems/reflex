[FILES]

Each line of the input file contains one waypoint, given as
space-separated floating point values.  The first value is the time to
be at that waypoint.  The second value is the blend time for that
waypoint.  The next four values are the orientation quaternion in xyzw
order.  The next three values are the translation xyz.

The output file contains a dense trajectory, with positions, and
optionally velocities at each timestep.  Each line contains the values
for a single timestep.  The first value is the time for that step.
The next four values are the orientation quaternion in xyzw order.
The next three values are the translation xyz.  If velocity is
requested, the next six values are the translational xyz then
rotational xyz velocities.
