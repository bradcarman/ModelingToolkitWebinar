using ModelingToolkit
using DifferentialEquations
using Plots

@parameters t
D = Differential(t)


# Connectors ----
# https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/connectors/connections/
@connector Port begin
    p(t)
    dm(t), [connect = Flow]
end

@connector Flange begin
    dx(t)
    f(t), [connect = Flow]
end


# Components ----
@mtkmodel Orifice begin
    @parameters begin
        Cₒ=2.7
        Aₒ=0.00094
        ρ₀=1000
    end
    @variables begin
        dm(t)=0
        p₁(t)
        p₂(t)
    end
    @components begin
        port₁ = Port()
        port₂ = Port()
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
    end
    @variables begin
        p(t)
        x(t)
        dm(t)=0
        f(t)
        dx(t)=0
        dr(t)=0
        (r(t) = dr ~ 0), [guess = 1000]
    end
    @components begin
        port = Port()
        flange = Flange()
    end
    @equations begin
        D(x) ~ dx
        D(r) ~ dr
        
        p ~ +port.p
        dm ~ +port.dm # mass is entering
        flange.f * direction  ~ -f 
        flange.dx * direction ~ dx 

        r ~ ρ₀*(1 + p/β)
        dm ~ (r*dx*A) + (dr*x*A)
        f ~ p * A
    end
end

@mtkmodel Mass begin
    @parameters begin
        m = 100
    end
    @variables begin
        dx(t)=0
        f(t)
    end
    @components begin
        flange = Flange()
    end
    @equations begin
        # connectors
        flange.dx ~ dx
        flange.f ~ -f
        
        # physics
        f ~ m*D(dx)
    end
end

@mtkmodel Actuator begin
    @parameters begin
        x=0.5
    end
    @components begin
        port₁ = Port()
        port₂ = Port()
        flange = Flange()

        vol₁ = Volume(;x, direction=-1)
        vol₂ = Volume(;x, direction=+1)
        mass = Mass()
        
    end
    @equations begin
        connect(port₁, vol₁.port)
        connect(port₂, vol₂.port)
        connect(vol₁.flange, vol₂.flange, mass.flange, flange)
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
        d = 1
    end
    @variables begin
        dx(t)
        f(t)
    end
    @components begin
        flange = Flange()
    end
    @equations begin
        # connectors
        flange.dx ~ dx
        flange.f ~ -f
        
        # physics
        f ~ d*dx
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

initsys = initializesystem(sys)





using Symbolics
eq = full_equations(sys)[7]
defs = ModelingToolkit.defaults(sys);
defs[sys.res₁.dm]
eq = substitute(eq, sys.res₁.dm => defs[sys.res₁.dm])
defs[sys.act.vol₁.r]
eq = simplify(substitute(eq, sys.act.vol₁.r => defs[sys.act.vol₁.r]))
defs[sys.act.vol₁.p′]
eq = substitute(eq, sys.act.vol₁.p′ => defs[sys.act.vol₁.p′])
defs[sys.act.p₁′]
eq = substitute(eq, sys.act.p₁′ => defs[sys.act.p₁′])
defs[sys.src.p′]
eq = substitute(eq, sys.src.p′ => defs[sys.src.p′])

eq = full_equations(sys)[7].rhs;
ModelingToolkit.fixpoint_sub(eq, defs)


using SimpleEuler: BackwardEuler
prob = ODEProblem(sys, [], (0,1e-12))

prob.f

sol = solve(prob, BackwardEuler(); dt=1e-12, adaptive=false)



for (s,v) in zip(states(sys), sol.u[2])
    println("$s => $v")
end

prob = ODEProblem(sys, [sys.act.vol₁.r => 1014.9], (0,0.1))
solve(prob; abstol=1, reltol=1)

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
