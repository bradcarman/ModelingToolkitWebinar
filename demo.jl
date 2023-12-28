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
    dm(t)=0, [connect = Flow]
end

@connector Flange begin
    dx(t)=0
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
        dm(t)=0
        p₁(t)=p′
        p₂(t)=p′
    end
    @components begin
        port₁ = Port(p=p′)
        port₂ = Port(p=p′)
    end
    begin
        u = dm/(ρ₀*Aₒ)
    end
    @equations begin
        dm ~ +port₁.dm
        dm ~ -port₂.dm
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
        dm(t)=0
        f(t)=p′ * A
        dx(t)=0
        (r(t) = dr ~ 0), [guess = 1000]
        dr(t)
    end
    @components begin
        port = Port(p=p′)
        flange = Flange(f=-p′ * A * direction)
    end
    @equations begin
        D(x) ~ dx
        D(r) ~ dr
        
        p ~ +port.p
        dm ~ +port.dm # mass is entering
        f ~ -flange.f * direction # force is leaving
        dx ~ flange.dx * direction

        r ~ ρ₀*(1 + p/β)
        dm ~ (r*dx*A) + (dr*x*A)
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
        dx(t)=0
        ẍ(t)=f′/m
    end
    @components begin
        flange = Flange(f=f′)
    end
    @equations begin
        D(x) ~ dx
        D(dx) ~ ẍ

        f ~ flange.f
        dx ~ flange.dx

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
        port₁ = Port(p=p₁′)
        port₂ = Port(p=p₂′)
        vol₁ = Volume(p′=p₁′, x′=x′,  direction=-1)
        vol₂ = Volume(p′=p₂′, x′=x′,  direction=+1)
        mass = Mass(f′=(p₂′ - p₁′)*A)
        flange = Flange(f=0)
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
        port = Port(p=p′)
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
        flange = Flange(f=0)
    end
    @equations begin
        flange.f ~ c*flange.dx
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

# @mtkbuild sys = System()
# prob = ODEProblem(sys, [sys.act.vol₁.r], (0,0.1))


@named odesys = System()
sys = structural_simplify(odesys)
sys_init = initializesystem(sys)


#=
function check_eqs(sys::ODESystem, ieq::Int)

    varmap = ModelingToolkit.defaults(sys)
    varmap = Dict(Symbolics.diff2term(ModelingToolkit.value(k)) => ModelingToolkit.value(varmap[k]) for k in keys(varmap))
    eqs = full_equations(sys)

    eq = eqs[ieq].rhs
    
    eq = ModelingToolkit.fixpoint_sub(eq, varmap)
    
    return eq
end


sub1 = varmap[sys.act.vol₁.r]
sub2 = varmap[sys.act.vol₁.p′]

sub1 = (varmap[sys.act.vol₁.ρ₀]*(1 + varmap[sys.act.p₁′]/varmap[sys.act.vol₁.β]))

-varmap[sys.src.p′] + (-varmap[sys.act.vol₁.β]*(varmap[sys.act.vol₁.ρ₀] - sub1)) / act₊vol₁₊ρ₀ + 0.5res₁₊Cₒ*res₁₊ρ₀*((res₁₊dm(t) / (res₁₊Aₒ*res₁₊ρ₀))^2)

eq = eqs[7].rhs

defs = ModelingToolkit.defaults(sys);
ModelingToolkit.fixpoint_sub(eq.arguments[1], defs)
ModelingToolkit.fixpoint_sub(eq.arguments[2], defs)
ModelingToolkit.fixpoint_sub(eq.arguments[3], defs)

@mtkbuild sys = System()

check_eqs(sys, 7)


unicodes = [
    "₊" => "_"
    "₁" => "1"
    "₂" => "2"
    "′" => "_p"
    "ₒ" => "o"
    "₀" => "0"
    "ρ" => "rho"
    "β" => "B"
]


function modelica_equations(sys::ODESystem)

    reps = [
        "Differential(t)" => "der"
        "~" => "="
        "(t)" => ""
    ]

    push!(reps, unicodes...)

    mo_eqs = String[]
    for eq in full_equations(sys)
        push!(mo_eqs, replace(string(eq), reps...) * ";")
    end

    return mo_eqs
end

function modelica_parameters(sys::ODESystem, prob::ODEProblem)

    
    mo_pars = String[]
    for (p,v) in zip(parameters(sys), prob.p)
        st = replace(string(p),unicodes...)
        push!(mo_pars, "parameter Real $st = $v;")
    end

    return mo_pars
end

function modelica_states(sys::ODESystem, prob::ODEProblem)

    reps = ["(t)"=>""]
    push!(reps, unicodes...)
    
    mo_sts = String[]
    for (p,v) in zip(states(sys), prob.u0)
        st = replace(string(p),reps...)
        push!(mo_sts, "Real $st(start = $v);")
    end

    return mo_sts
end

mo_eqs = modelica_equations(sys)
foreach(println, mo_eqs)

mo_pars = modelica_parameters(sys, prob)
foreach(println, mo_pars)

mo_sts = modelica_states(sys, prob)
foreach(println, mo_sts)

prob = ODEProblem(sys, [], (0, 0.1), [])
=#

# https://docs.sciml.ai/DiffEqDocs/stable/solvers/dae_solve/#Initialization-Schemes
sol = solve(prob, ImplicitEuler(nlsolve = NLNewton(check_div=false, always_new=true)))
sol = solve(prob, Rodas4(); reltol=1e-1, abstol=1e-1)
# sol′ = solve(prob, Rodas5P(), reltol=1e-8, abstol=1e-8, initializealg = ShampineCollocationInit())

# velocity comparison (incompressible vs. compressible)
plot(sol, idxs=[sys.act.mass.dx]; ylabel="velocity [m/s]")
plot!(sol_ic, idxs=[ẋ])


# What's Next --> Using the ModelingToolkitStandardLibrary
# https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/
# RC Circuit
# https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/tutorials/rc_circuit/
# DC Motor
# https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/tutorials/dc_motor_pi/
