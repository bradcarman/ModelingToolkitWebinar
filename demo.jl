using ModelingToolkit
using DifferentialEquations
using Plots

# ------------------------------------------------
# Part 1: Steady State Modeling ------------------
# ------------------------------------------------
vars = @variables p₁=300e5 p₂=0e5 Aₒ=0.001
pars = @parameters A=0.1 ẋ=1 c=1000 pₛ=300e5 pᵣ=0 ρ=1000 Cₒ=2.7 m=100 ẍ=0

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
    sol = solve(prob′)
    push!(orifices, sol[Aₒ])
end
plot(velocity_limits, orifices; xlabel="velocity limit [m/s]", ylabel="orifice size [m^2]")



# ------------------------------------------------
# Part 2: Dynamic Modeling (DAEs) ----------------
# ------------------------------------------------
@parameters t
D = Differential(t)

pars = @parameters A=0.1 pₛ=300e5 pᵣ=0 ρ=1000 C₀=2.7 m=100 Aₒ=0.00094 c=1000
vars = @variables x(t)=0 ẋ(t)=0 p₁(t)=300e5 p₂(t)=0e5 ẍ(t)=(p₂-p₁)*A

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

# Connectors ----
# https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/connectors/connections/
@connector Port begin
    p(t)
    ṁ(t)=0, [connect = Flow]
end

@connector Flange begin
    ẋ(t)=0
    f(t), [connect = Flow]
end


# Components ----
@mtkmodel Orifice begin
    @parameters begin
        Cₒ=2.7
        Aₒ=0.00094
        ρ₀=1000
        p′=0
    end
    @variables begin
        ṁ(t)=0
        p₁(t)=p′
        p₂(t)=p′
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
        
        p₁ - p₂ ~ (1/2)*ρ₀*u^2*Cₒ
    end
end

@mtkmodel Volume begin
    @parameters begin
        A=0.1
        ρ₀=1000
        β=2e9
        direction=+1
        p′
        x′
    end
    @variables begin
        p(t)=p′
        x(t)=x′
        ṁ(t)=0
        f(t)=p′ * A
        ẋ(t)=0
        r(t)=ρ₀*(1 + p′/β)
        ṙ(t)=0
    end
    @components begin
        port = Port()
        flange = Flange()
    end
    @equations begin
        D(x) ~ ẋ
        D(r) ~ ṙ
        
        p ~ +port.p
        ṁ ~ +port.ṁ # mass is entering
        f ~ -flange.f * direction # force is leaving
        ẋ ~ flange.ẋ * direction

        r ~ ρ₀*(1 + p/β)
        ṁ ~ (r*ẋ*A) + (ṙ*x*A)
        f ~ p * A
    end
end

@mtkmodel Mass begin
    @parameters begin
        m = 100
        f′
    end
    @variables begin
        f(t)=f′
        x(t)=0
        ẋ(t)=0
        ẍ(t)=f′/m
    end
    @components begin
        flange = Flange()
    end
    @equations begin
        D(x) ~ ẋ
        D(ẋ) ~ ẍ

        f ~ flange.f
        ẋ ~ flange.ẋ

        m*ẍ ~ f
    end
end

@mtkmodel Actuator begin
    @parameters begin
        p₁′
        p₂′
    end
    begin #constants
        x′=0.5
        A=0.1
    end
    @components begin
        port₁ = Port()
        port₂ = Port()
        vol₁ = Volume(p′=p₁′, x′=x′,  direction=-1)
        vol₂ = Volume(p′=p₂′, x′=x′,  direction=+1)
        mass = Mass(f′=(p₂′ - p₁′)*A)
        flange = Flange()
    end
    @equations begin
        connect(port₁, vol₁.port)
        connect(port₂, vol₂.port)
        connect(vol₁.flange, vol₂.flange, mass.flange, flange)
    end
end

@mtkmodel Source begin
    @parameters begin
        p′
    end
    @components begin
        port = Port()
    end    
    @equations begin
        port.p ~ p′
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
        flange.f ~ c*flange.ẋ
    end
end


@mtkmodel System begin
    @components begin
        res₁ = Orifice(p′=300e5)
        res₂ = Orifice(p′=0)
        act = Actuator(p₁′=300e5, p₂′=0)
        src = Source(p′=300e5)
        snk = Source(p′=0)
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
prob = ODEProblem(sys, [], (0, 0.1), [])

# https://docs.sciml.ai/DiffEqDocs/stable/solvers/dae_solve/#Initialization-Schemes
NEWTON = NLNewton(check_div = false, always_new = true, max_iter = 1000, relax = 4 // 10)
sol = solve(prob, ImplicitEuler(nlsolve = NEWTON); initializealg = NoInit(), dt = 1e-6, adaptive = false)
sol′ = solve(prob, ImplicitEuler(), reltol=1e-8, abstol=1e-8, initializealg = ShampineCollocationInit())

println("NoInit() adaptive=false")
println("state = given => used")
for (s, u0, st) in zip(sol[1], prob.u0, states(sys))
    println("$st = $u0 => $s")
end

println("ShampineCollocationInit() adaptive=true")
println("state = given => used")
for (s, u0, st) in zip(sol′[1], prob.u0, states(sys))
    println("$st = $u0 => $s")
end


plot(sol, idxs=[sys.act.vol₁.p, sys.act.vol₂.p])
plot(sol, idxs=[sys.act.vol₁.x, sys.act.vol₂.x])

plot(sol, idxs=[sys.act.vol₁.ṁ , sys.act.vol₂.ṁ ])

plot(sol, idxs=[sys.act.mass.x])
plot!(sol_ic, idxs=[x])

plot(sol, idxs=[sys.act.mass.ẋ])
plot!(sol_ic, idxs=[ẋ])

plot(sol, idxs=[sys.act.mass.ẍ])
plot!(sol_ic, idxs=[ẍ])


# What's Next --> Using the ModelingToolkitStandardLibrary
# https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/