include("hex_to_orth.jl")
include("utilities/particle_moments.jl")
using Test
p=Parameters(eps_c=0.01,anisotropy=3.0); grid=make_grid(16,16); B=elastic_kernels(grid,p)
c,eta=initial_fields(16)
eta = [e .+ 0.15a for (a,e) in enumerate(eta)]
dc,de=bulk_derivatives(c,eta,p); ec,ee=elastic_derivatives(c,eta,B)
vc=randn(MersenneTwister(4),size(c)); ve=[randn(MersenneTwister(a),size(c)) for a in 1:3]
eps=1e-6
numerical=(total_energy(c.+eps.*vc,[eta[a].+eps.*ve[a] for a in 1:3],grid,B,p)-
           total_energy(c.-eps.*vc,[eta[a].-eps.*ve[a] for a in 1:3],grid,B,p))/(2eps)
mu=dc.+ec.+real.(ifft(2p.kappa_c.*grid.k2.*fft(c)))
deta=[de[a].+ee[a].+real.(ifft(2p.kappa_eta.*grid.k2.*fft(eta[a]))) for a in 1:3]
analytical=(sum(mu.*vc)+sum(sum(deta[a].*ve[a]) for a in 1:3))*grid.dx*grid.dy
@test isapprox(numerical,analytical;rtol=1e-6,atol=1e-6)
for i in axes(B,1),j in axes(B,2)
    @test minimum(eigvals(Symmetric(B[i,j,:,:]))) > -1e-8
end
_,_,history=run_example()
@test maximum(abs.(history[:,2].-history[1,2])) < 1e-12
@test maximum(diff(history[:,3])) < 1e-10
println("Energy derivative, kernel positivity, finite fields and mass checks passed.")
println("Mean drift: ",maximum(abs.(history[:,2].-history[1,2])))
println("Energy: ",history[1,3]," → ",history[end,3])
x=collect(-50:50)
circle=[i^2+j^2 <= 20^2 ? 1.0 : 0.0 for i in x,j in x]
r=particle_moments(circle,circle)
@test r.shape_parameter < 1e-12
@test isnan(r.long_axis_angle_degrees)
ellipse=[(i/30)^2+(j/15)^2 <= 1 ? 1.0 : 0.0 for i in x,j in x]
r=particle_moments(ellipse,ellipse)
@test abs(r.shape_parameter-0.6) < 0.02
shifted=circshift(ellipse,(5,8))
@test isapprox(particle_moments(shifted,shifted).principal_moments,r.principal_moments)
@test isapprox(particle_moments(ellipse,ellipse;dx=2,dy=2).principal_moments,16 .*r.principal_moments)
@test_throws ErrorException particle_moments(zeros(4,4),zeros(4,4))
println("Circle, ellipse, translation, length scaling and empty-mask utility checks passed.")
println("Ellipse shape parameter: ",r.shape_parameter)

# Known isotropic, in-plane dilatational kernel checks the elastic construction.
isotropic = Parameters(eps_c=0.01, eps_eta=0.0, anisotropy=1.0)
Bi = elastic_kernels(grid,isotropic)
C = elastic_tensor(isotropic)
lambda = C[1,1,2,2]
shear = C[1,2,1,2]
expected = 4shear*(lambda+shear)/(lambda+2shear)*isotropic.eps_c^2
for i in axes(Bi,1), j in axes(Bi,2)
    i == 1 && j == 1 && continue
    @test isapprox(Bi[i,j,1,1],expected;rtol=1e-12)
end
println("Known isotropic dilatational kernel check passed.")
