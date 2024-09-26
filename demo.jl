using ModelingToolkit
using DifferentialEquations
using Plots

# ------------------------------------------------
# Part 1: Steady State Modeling ------------------
# ------------------------------------------------
pars = @parameters A=0.1 ẋ=1 c=1000 pₛ=300e5 pᵣ=0 ρ=1000 Cₒ=2.7 m=100 ẍ=0
vars = @variables p₁=300e5 p₂=0e5 Aₒ=0.001

# symbolic expressions
u = ẋ * (A/Aₒ)

# equations
eqs = [
    pₛ - p₁ ~ (1/2)*ρ*u^2*Cₒ
    p₂ - pᵣ ~ (1/2)*ρ*u^2*Cₒ
    m*ẍ ~ (p₂ - p₁)*A - c*ẋ
]

@named nlsys = NonlinearSystem(eqs, vars, pars)
sys = structural_simplify(nlsys)
prob = NonlinearProblem(sys, [], []) # [initial conditions], [parameters] 
sol = solve(prob)

sol[Aₒ] #<-- solution!

# how to quickly make a new solution 
orifices = []
velocity_limits = 1.0:0.1:2.0
for velocity_limit in velocity_limits
    prob′ = remake(prob; p=[ẋ => velocity_limit])
    sol′ = solve(prob′)
    push!(orifices, sol′[Aₒ])
end
plot(velocity_limits, orifices; xlabel="velocity limit [m/s]", ylabel="orifice size [m^2]")





# ------------------------------------------------
# Part 2: Dynamic Modeling (DAEs) ----------------
# ------------------------------------------------
@parameters t
D = Differential(t)

pars = @parameters A=0.1 pₛ=300e5 pᵣ=0 ρ=1000 C₀=2.7 m=100 Aₒ=0.00094 c=1000
vars = @variables x(t)=0 ẋ(t)=0 p₁(t)=300e5 p₂(t)=0e5 ẍ(t)=(p₂-p₁)*A/m

# symbolic expressions
u = ẋ * (A/Aₒ)

# equations
eqs = [
    D(x) ~ ẋ
    D(ẋ) ~ ẍ

    pₛ - p₁ ~ (1/2)*ρ*u^2*C₀
    p₂ - pᵣ ~ (1/2)*ρ*u^2*C₀

    m*ẍ ~ (p₂-p₁)*A - c*ẋ
]

@named odesys = ODESystem(eqs, t, vars, pars)
sys = structural_simplify(odesys)
prob = ODEProblem(sys, [], (0.0, 0.0001), [])
sol = solve(prob)

# explain sol object...
plot(sol.t, sol[x]; marker=:circle, ylabel="position [m]")
plot(sol, idxs=[x]; ylabel="position [m]")
plot(sol, idxs=[ẋ]; ylabel="velocity [m/s]")
plot(sol, idxs=[ẍ]; ylabel="acceleration [m/s^2]")
plot(sol, idxs=[p₁, p₂]; ylabel="pressure [Pa]")

# for comparison with compressible system
prob′ = remake(prob, tspan=(0, 0.1))
sol_ic = solve(prob′)




# ------------------------------------------------
# Part 3: Component Based Modeling ---------------
# ------------------------------------------------
regPow(x, a, delta = 0.01) = x * (x * x + delta * delta)^((a - 1) / 2);
regRoot(x, delta = 0.01) = regPow(x, 0.5, delta)

# Connectors ----
# https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/connectors/connections/
@connector Port begin
    p(t), [guess=0]
    ṁ(t), [guess=0, connect = Flow]
end

@connector Flange begin
    x(t), [guess=0]
    f(t), [guess=0, connect = Flow]
end


# Components ----
@mtkmodel Orifice begin
    @parameters begin
        Cₒ=2.7
        Aₒ=0.00094
        ρ₀=1000
    end
    @variables begin
        ṁ(t), [guess=0]
        p₁(t), [guess=1]
        p₂(t), [guess=1]
    end
    @components begin
        port₁ = Port()
        port₂ = Port()
    end
    begin
        u = ṁ/(ρ₀*Aₒ)
    end
    @equations begin
        ṁ ~ +port₁.ṁ
        ṁ ~ -port₂.ṁ
        p₁ ~ port₁.p
        p₂ ~ port₂.p
        
        # p₁ - p₂ ~ (1/2)*ρ₀*u^2*Cₒ

        u ~ regRoot( 2*(p₁ - p₂)/(ρ₀*Cₒ) )

    end
end

@mtkmodel Volume begin
    @parameters begin
        A=0.1
        ρ₀=1000
        β=2e9
        direction=+1
        L=0.5
    end
    @variables begin
        p(t), [guess=1]
        x(t), [guess=0]
        m(t), [guess=1]
        ṁ(t), [guess=0]
        f(t), [guess=1]
        ẋ(t), [guess=0]
        r(t), [guess=1]
        ṙ(t), [guess=0]
    end
    @components begin
        port = Port()
        flange = Flange()
    end
    @equations begin
        D(x) ~ ẋ
        D(r) ~ ṙ
        D(m) ~ ṁ
        
        p ~ +port.p
        ṁ ~ +port.ṁ # mass is entering
        f ~ -flange.f * direction # force is leaving
        ẋ ~ +D(flange.x) * direction

        r ~ ρ₀*(1 + p/β)
        m ~ r*(x+L)*A
        f ~ p * A
    end
end

@mtkmodel Mass begin
    @parameters begin
        m = 100
    end
    @variables begin
        f(t), [guess=1]
        x(t), [guess=0]
        ẋ(t), [guess=0]
        ẍ(t), [guess=0]
    end
    @components begin
        flange = Flange()
    end
    @equations begin
        D(x) ~ ẋ
        D(ẋ) ~ ẍ

        f ~ flange.f
        x ~ flange.x

        m*ẍ ~ f
    end
end

@mtkmodel Actuator begin
    @parameters begin
        A=0.1
    end
    @variables begin
        x(t), [guess=0]
    end
    @components begin
        port₁ = Port()
        port₂ = Port()
        vol₁ = Volume(;A,  direction=-1)
        vol₂ = Volume(;A,  direction=+1)
        mass = Mass()
        flange = Flange()
    end
    @equations begin
        connect(port₁, vol₁.port)
        connect(port₂, vol₂.port)
        connect(vol₁.flange, vol₂.flange, mass.flange, flange)

        x ~ mass.x
    end
end

@mtkmodel Source begin
    @parameters begin
        p
    end
    @components begin
        port = Port()
    end    
    @equations begin
        port.p ~ p
    end
end

@mtkmodel Damper begin
    @parameters begin
        c = 1000
    end
    @components begin
        flange = Flange()
    end
    @equations begin
        flange.f ~ c*D(flange.x)
    end
end


@mtkmodel System begin
    @components begin
        res₁ = Orifice()
        res₂ = Orifice()
        act = Actuator()
        src = Source(p=300e5)
        snk = Source(p=0)
        dmp = Damper()
    end
    @equations begin
        connect(src.port, res₁.port₁)
        connect(res₁.port₂, act.port₁)
        connect(act.port₂, res₂.port₁)
        connect(res₂.port₂, snk.port)
        connect(dmp.flange, act.flange)
    end
end

@mtkbuild sys = System()

ρ₀=1000
β=2e9
p1=300e5
p2=0
u0 = unknowns(sys) .=> [
    0
    0
    0.5*0.1*ρ₀*(1 + p1/β)
    0
    0.5*0.1*ρ₀*(1 + p2/β)
    0
    0
    0
    ρ₀*(1 + p1/β)
    ρ₀*(1 + p2/β)
]

include("convert_to_modelica.jl")
convert_to_modelica(sys, Dict(u0))



initialization_eqs = [
    sys.act.x ~ 0
    D(sys.act.x) ~ 0
    
    sys.act.vol₁.x ~ 0
    D(sys.act.vol₁.m) ~ 0
    
    sys.act.vol₂.x ~ 0
    D(sys.act.vol₂.m) ~ 0
]



initsys = ModelingToolkit.generate_initializesystem(sys; initialization_eqs)
structural_simplify(initsys)

prob = ODEProblem(sys, u0, (0, 0.1), []; initialization_eqs)
sol = solve(prob)


# velocity comparison (incompressible vs. compressible)
plot(sol, idxs=[sys.act.mass.ẋ]; ylabel="velocity [m/s]", label="Compressible")
plot!(sol_ic, idxs=[ẋ], label="Incompressible")


# What's Next --> Using the ModelingToolkitStandardLibrary
# https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/
# RC Circuit
# https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/tutorials/rc_circuit/
# DC Motor
# https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/tutorials/dc_motor_pi/
