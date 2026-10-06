using LinearAlgebra

# I store x along rows and y along columns, and measure one isolated particle.
# I first shift a boundary-crossing particle into the box before measuring it.
function particle_moments(c, eta; dx=1.0, dy=1.0, threshold=0.5)
    size(c) == size(eta) || error("Field sizes differ.")
    dx > 0 && dy > 0 || error("Grid spacing must be positive.")
    all(isfinite,c) && all(isfinite,eta) || error("Fields must be finite.")
    points = findall((c .>= threshold) .& (eta .>= threshold))
    length(points) > 1 || error("Select a particle containing at least two points.")
    x = [(p[1]-1)*dx for p in points]
    y = [(p[2]-1)*dy for p in points]
    center = (sum(x)/length(x), sum(y)/length(y))
    x .-= center[1]; y .-= center[2]
    Ixx = sum(y.^2)*dx*dy
    Iyy = sum(x.^2)*dx*dy
    Ixy = -sum(x.*y)*dx*dy
    tensor = [Ixx Ixy; Ixy Iyy]
    principal = eigen(Symmetric(tensor))
    small, large = principal.values
    large+small > 0 || error("Zero second moment.")
    shape_parameter = (large-small)/(large+small)
    axis = principal.vectors[:,1]
    angle = shape_parameter < 1e-12 ? NaN : mod(atand(axis[2],axis[1]),180)
    area = length(points)*dx*dy
    return (; center, area, equivalent_radius=sqrt(area/pi), tensor,
        principal_moments=principal.values, shape_parameter, long_axis_angle_degrees=angle)
end

# I read the original interleaved real and imaginary doubles, with y changing fastest.
function read_c_field(path,nx,ny)
    values = reinterpret(Float64,read(path))
    length(values) == 2nx*ny || error("File size does not match nx and ny.")
    return permutedims(reshape(values[1:2:end],ny,nx))
end
