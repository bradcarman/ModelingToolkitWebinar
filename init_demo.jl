using ModelingToolkit

@parameters t
D = Differential(t)

@mtkmodel Flange begin
    @variables begin
        dx(t), [guess = 0]
        f(t), [guess = 0, connect=Flow]
    end
end

@mtkmodel Mass begin
    @parameters begin
        m = 100
    end
    @variables begin
        dx(t)
        f(t)=0
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

@mtkmodel Damper begin
    @parameters begin
        d = 1
    end
    @variables begin
        dx(t), [guess = 0]
        f(t), [guess = 0]
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
        mass = Mass(;dx=100,m=10)
        damper = Damper(;d=1)
    end
    @equations begin
        # connect(mass.flange, damper.flange)
        mass.flange.dx ~ damper.flange.dx
        0 ~ mass.flange.f + damper.flange.f
    end
end

@named odesys = System()
equations(odesys)

@mtkbuild sys = System()
# prob = ODEProblem(sys, [], (0,1))
# sol = solve(prob)

isys = initializesystem(sys)
iprob = NonlinearProblem(isys, []; check_length=false)
eqs = equations(isys)

popat!(eqs, 4)

@named isys = NonlinearSystem(eqs, states(isys), [])
iprob = NonlinearProblem(isys, iprob.u0)
isol = solve(iprob)