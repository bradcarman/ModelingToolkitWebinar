using ModelingToolkit
using DifferentialEquations
using Plots

@parameters t
D = Differential(t)


# Connectors ----
# https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/connectors/connections/
@connector Port begin
    p(t)
    dm(t), [connect=Flow, guess=0]
end

@connector Flange begin
    dx(t), [guess=0]
    f(t), [connect=Flow, guess=0]
end


# Components ----
@mtkmodel Orifice begin
    @parameters begin
        Cₒ=2.7
        Aₒ=0.00094
        ρ₀=1000
        p
    end
    @variables begin
        dm(t)=0
        p₁(t)=p
        p₂(t)=p
    end
    @components begin
        port₁ = Port(;p=p₁)
        port₂ = Port(;p=p₂)
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
        f(t), [guess=0]
        dx(t)=0
        dr(t)=0
        (r(t) = dr ~ 0), [guess = 1000]
    end
    @components begin
        port = Port(;p)
        flange = Flange(;f)
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
        f(t), [guess=0]
    end
    @components begin
        flange = Flange(;f)
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
    @variables begin
        p₁(t)
        p₂(t)
    end
    @components begin
        port₁ = Port(;p=p₁)
        port₂ = Port(;p=p₂)
        flange = Flange()

        vol₁ = Volume(;x, p=p₁, direction=-1)
        vol₂ = Volume(;x, p=p₂, direction=+1)
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
        port = Port(;p)
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
        dx(t)=0
        f(t), [guess=0]
    end
    @components begin
        flange = Flange(;f)
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
    @parameters begin
        p₁=300e5
        p₂=0
    end
    @components begin
        res₁ = Orifice(;p=p₁)
        res₂ = Orifice(;p=p₂)
        act = Actuator(;p₁, p₂)
        src = Source(p=300e5)
        snk = Source(p=0)
        damper = Damper()
    end
    @equations begin
        connect(src.port, res₁.port₁)
        connect(res₁.port₂, act.port₁)
        connect(act.port₂, res₂.port₁)
        connect(res₂.port₂, snk.port)
        connect(damper.flange, act.flange)
    end
end

@mtkbuild sys = System()

initsys = initializesystem(sys)
