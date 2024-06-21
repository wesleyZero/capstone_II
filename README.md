# Capstone II
**UC Santa Barbara Chemical Engineering 2024**

## Table of Contents 

 * [Description of the Project](#Description)
 * [Overview of the project](https://github.com/wesleyZero/capstone_II?tab=readme-ov-file#overview-of-the-project)
 * [How do you simulate a chemical process?](https://github.com/wesleyZero/capstone_II?tab=readme-ov-file#so-how-do-you-simulate-a-chemical-reactor)
    * [Process Flow Diagram](https://github.com/wesleyZero/capstone_II?tab=readme-ov-file#process-flow-diagram)
    * [The Primary functions (higher level logoic)](https://github.com/wesleyZero/capstone_II?tab=readme-ov-file#the-primary-functions-the-higher-level-logic)
    * [Reaction Chemistry](https://github.com/wesleyZero/capstone_II?tab=readme-ov-file#reaction-chemistry)
    * [Chemical Kinetics](https://github.com/wesleyZero/capstone_II?tab=readme-ov-file#chemical-kinetics)
       * [Rate Constant](https://github.com/wesleyZero/capstone_II?tab=readme-ov-file#rate-constant)
       * [Reaction Rate](https://github.com/wesleyZero/capstone_II?tab=readme-ov-file#reaction-rate)
    * [Solving the system of equations](https://github.com/wesleyZero/capstone_II?tab=readme-ov-file#solving-the-system-of-non-linear-equations)
    * [The recycle stream](https://github.com/wesleyZero/capstone_II?tab=readme-ov-file#the-recycle-stream)
    * [Generate Plots to Answer the Original Questions](https://github.com/wesleyZero/capstone_II?tab=readme-ov-file#genenate-plots-to-answer-the-original-questions)
 * [Nomenclature](#Nomenclature)
    * [Supercritical Fluids](#Supercritical_fluids)
 * [Links](#Links)

## Description
This description and readme is a _massive_ simplification of our chemical engineering capstone for the intended audience of software engineers, with the entire emphasis being of how the code was written. For the full detailed report (with the intended audience of chemical engineers) please see the [PDF of our final Report](https://github.com/wesleyZero/capstone_II/blob/main/readme/pdf/capstone2_DMC.pdf). We (my capstone group) performed a detailed techno-economic analysis and design for a proposed 100 kiloton per year (kta) plant producing polymer-grade Dimethyl Carbonate (DMC) for optical applications. Our economic analysis included (but is certainly not limited to) calculations of net present value (NPV), internal rate of return (IRR), total capital investment (TCI), energy consumption, carbon emissions, hazard ananylsis, and _much_ more. Those details will be ommitted in this readme since the focus of the following will be the code. 

# Overview of the project

**The end goal of this simulation is to answer just a few questions.** 
- What size of a reactor are we going to use?
- Are we going to operate the reactor isothermally? or isobarically? (i.e. at constant temp or pressure?)
     - If isothermal, what pressure are we going to operate at? (note: we were only allowed one isothermal temperature)
     - If isobaric, what temperature are we going to operate at? (note: we were only allowed one isobaric pressure to operate at)
- What is the economic value of the reactor that we designed? how do we choose conditions to optimize the economics of the process?

**How will we answer these questions?**
- Look at the reaction chemistry (what reactions can occur?, what species are involved?)
- Evaluate the chemical kinetics (given a set of conditions, how much of each reaction will occur?)
  
Whats important to know is **the end goal is to produce a set of graphs**. These graphs are going be be **functions of the residence time** (the mean time a chemical species spends in the reactor, which is proportional to the reactor volume, symbol: tau), **conversion** (how much of the limiting species gets converted into products, symbol: chi), **temperature**, and **pressure**. The other important thing to note, is that we initially didn't know if an isothermal (constant temp) or isobaric (constant pressure) model was going to work. So we had to try out both, you will see this in the code.

# So, How do you simulate a chemical reactor?

You start off with a process flow diagram, this process flow diagram shows all processes that are included in the chemical processing plant that we are designing. It's good to have an idea of what the overall picture is. 

## Process Flow Diagram 
![Process Flow Diagram](https://github.com/wesleyZero/capstone_II/blob/main/readme/img/process_flow_diagram.jpeg)

**You can see above** all of the unit opererations (reactors, separators, distillation towers, etc.). **This readme will only be focusing on the reactor** near the top center-left of the diagram. **Why?** Well, in short chemical separations involve lots of "fudge factors" i.e. correlation factors and empirical models that account for the differences between simulation and reality. With a separation process this complex, it's not very useful to simulate the entire separation since the amount of error that will accumulate throughout the process is so large that it's neccessary to employ professional process simulation software like the one we used (Aspen HYSYS). Even the process simulation software has a large degree of error and requires a very educated chemical engineer to use correctly, and _even then_ there will be error.

<p align="center">
  <img src="https://github.com/wesleyZero/capstone_II/blob/main/readme/img/reactor.png" alt="CSTR Reactor">
</p>

The reactor we are using for this process is called a CSTR (constantly stirred tank reactor). It's kinda like if you imagine a witch brewing an evil potion or something like that 🧙‍♂️, however there are pipes that are flowing in and out of the brew constantly flowing reactants in and products out, there are lots of differential equations (which turn algebraic at steady state), and in this particular case some of the reactants are super critical fluids! Simple right?! **But Wes...What the heck is a super-critical fluid?** Read here <a href="https://github.com/wesleyZero/capstone_II/blob/main/README.md#supercritical_fluidsf" target="_blank">Supercritical Fluids</a> to learn what this is, if you wanna know! If you don't care what a supercritical fluid is, just know that it has variable density and that changes the concentration and reaction rate of everything in the reactor

## The primary Functions (the higher level logic)
[level3_flowrates](https://github.com/wesleyZero/capstone_II/blob/main/run_dmc.m#L198)
```matlab
function [F_fresh, F_rxtr, F_out, R, V_rxtr] = level3_flowrates(tau, temp, P, opt)
    F_fresh = NaN; F_rxtr = NaN; F_out = NaN; R = NaN;
    user = get_user_inputs(); 
    flow_fxns = flowrate_fxns();
    rxtr_fxns = reactor_fxns();

    F_basis = flow_fxns.get_basis_feed_flowrates();
    [F_fresh, F_rxtr, F_out, R, V_rxtr] = rxtr_fxns.get_reactor_flows(F_basis, temp, P, opt, tau);
    

end
```

This is the core function (the root of the call stack) as far as you are concerned to keep things simple. What we are going to do is call this function for a bunch of different **residence times (tau)** and conditions (T, P, isobaric/isothermal). We will plot the outputs of this and evaluate the economics of this functions outputs to create the graphs and find answers to the questions we want to answer. **Important note about residence time:** For a given flowrate, the residence time is proportional to the volume of the reactor. We create functions of residence time though because it's an intrinsic property of the system, and we want our calculations to be arbitrarily scalable. 

[get_reactor_flows](https://github.com/wesleyZero/capstone_II/blob/main/reactor_fxns.m#L9)
```matlab
function [F_fresh, F_real_feed, F_real_eff, R, V_rxtr] = get_reactor_flows(F_real_feed_basis, T, P, opt, tau)
    F_fresh = NaN; F_real_feed = NaN; F_real_eff = NaN; R = NaN;
    % basis calculations for the real reactor 
    q_tot.basis = get_total_volumetric_flowrate(F_real_feed_basis, T, P, opt);
    V_rxtr.basis = q_tot.basis * tau;
    C_out = get_reactor_effluent_concentrations(F_real_feed_basis, T, P, opt, tau);

    if any(imag(C_out) ~= 0)
        % disp('ERROR : Complex valued concentrations');
        return 
    elseif any(real(C_out) < 0)
        % disp('ERROR : Negative valued concentrations');
        return
    else
        % disp('Valid solution!!!!!!!!!!')
        C_out = get_concentration_struct(C_out);
    end

    F_real_eff_basis = conc_to_flowrate(C_out, q_tot.basis);
        % assumption : liquid flow has no change in vol in effluent 

    % Plant Scale Calculations 
    [F_fresh, F_real_feed, F_real_eff, R] = get_plant_flowrates(F_real_feed_basis, F_real_eff_basis);
    scale_factor = get_scale_factor(F_real_eff_basis);
    V_rxtr.plant = V_rxtr.basis * scale_factor ; 

    V_rxtr = V_rxtr.plant;

end
```

This function (called by the one before) is also a core part of the logic. We start out with a **basis** for our calculations, very analogous to a vector that can be scaled. Our basis contains the intrinsic information about our feed composition. It's of an abiturary size, **since we will scale it at the end to satisfy the requirements of how much product we need to generate.** We turn the molar flowrate into a volumetric flowrate, using functions of density. This will change due to the supercritical fluid density changing. Then we can calculate the volume of the reactor and get the effluent (output) concentrations, by knowing how long reactions were allowed to take place (tau) and what the concentrations of all the species where at the provided temperature. We then make sure we have valid solution, scale it to make sure we are outputting 100 kta of DMC, and return the values we get. 


## Reaction Chemistry

First, we need to know what chemical species are reacting, and how much. So we look at the reaction chemistry. 

![](https://github.com/wesleyZero/capstone_II/blob/main/readme/img/reaction_chemistry.png)

## Chemical Kinetics
Chemical kinetics, simply means "how fast is this reaction happening?" i.e. how much will the concentration of a species change per second?. The equations below are how we mathematically determine the reaction rate. The first reaction is assumed to be instantaneous. To determine this, we must first look at the rate constant. 

<p align="center">
  <img src="https://github.com/wesleyZero/capstone_II/blob/main/readme/img/kinetics_2.png">
</p>

<p align="center">
  <img src="https://github.com/wesleyZero/capstone_II/blob/main/readme/img/kinetics_3.png">
</p>

### Rate constant

The **rate constant** is a constant that tells you how fast a particular reaction will happen, for a given concentration of reactants. The rate constant is a function of temperature, because at higher tempertures we get more chemical collisions and a higher concentration of species with enough energy to react upon a collision.

[get_isobaric_rate_constant](https://github.com/wesleyZero/capstone_II/blob/main/reactor_fxns.m#L323)
```matlab
function k = get_isobaric_rate_constant(reaction, T)
    % input: 
    %   T [ C ]
    % output:
    %   k [mol / L s ]

    const = get_constants();
    thermo = const.thermo;
    T = const.units.temperature.c_to_k(T); % [ K ]

    switch reaction
        case '2f'
            k = 6.69 * 10^2 * exp(-37200 / (thermo.R * T));
        case '2r'
            k = 1.19 * 10^4 * exp(-53700 / (thermo.R * T));
        case '3'
            k = 1.89 * 10^6 * exp(-82400 / (thermo.R * T));
        otherwise
            k = NaN;
            disp("ERROR: get_isobaric_rate_constant(): invalid reaction option");
    end
end
```

[get_isothermal_rate_constant](https://github.com/wesleyZero/capstone_II/blob/main/reactor_fxns.m#L346)
```matlab
function k = get_isothermal_rate_constant(reaction, T, P)
    % input:
    %   P [ bar ]
    % These functions are from the research paper
    rho = get_supercritical_c02_density(T, P, 'isothermal');
    switch reaction
        case '2f'
            if rho > 246.82 % [g / L]
                k = (2.486 * 10^(-2)) - (4.943 * (10^(-5)) * rho);
            else
                k = (1.362 * 10^(-2)) - (1.569 * (10^(-6)) * rho);
            end
        case '2r'
            k = 0.01486 * rho^(-0.873);
        case '3'
            k = 3.014 * (10^(-4)) * exp(-5.99 * (10^(-3)) * rho);
        otherwise
            k = NaN;
            disp("ERROR: get_isothermal_rate_constant(): invalid reaction option")
    end
end
```

### Reaction Rate

Once we have the rate constants, we can calculate the reaction rate. 

[get_reaction_rate](https://github.com/wesleyZero/capstone_II/blob/main/reactor_fxns.m#L207)
```matlab
function r = get_reaction_rate(C, reaction, T, P, opt)
    % input:
    %   opt = 'isothermal' or 'isobaric'

    k = get_rate_constant(reaction, T, P, opt);
    switch reaction
        case '2f'
            r = k * (C.ethylene_carbonate)^0.8;
        case '2r'
            r = k * C.dimethyl_carbonate * C.ethylene_carbonate;
        case '3'
            r = k * C.ethylene_carbonate;
        otherwise
            r = NaN;
            disp("ERROR: get_reaction_rate(): invalid reaction option")
    end
end
```

### Solving the system of non-linear equations

Now that we have the reaction rates we can use those reaction rates to solve a system of non-linear equations. This solution to this system will give us the concentrations coming out of the reactor. 

[sys_of_eqns](https://github.com/wesleyZero/capstone_II/blob/main/reactor_fxns.m#L177)
```matlab
function eqn = sys_of_eqns(C, params)
    T = params.T;
    P = params.P;
    opt = params.opt;
    Ci0 = params.Ci0;
    tau = params.tau;
    C_struct = get_concentration_struct(C);
    r = get_all_reaction_rates(C_struct, T, P, opt);

    r.ec =  r.r2r - r.r2f - r.r3; 
    r.meoh = (2 * r.r2r) - (2 * r.r2f) - r.r3;
    r.co2 = r.r3;
    r.dmc = r.r2f - r.r2r;
    r.eg = r.r2f - r.r2r; 
    r.me = r.r3; 

    eqn(1) = Ci0.ethylene_carbonate - C(1) + (tau * r.ec);
    eqn(2) = Ci0.methanol - C(3) + (tau * r.meoh);
    eqn(3) = Ci0.carbon_dioxide - C(4) + (tau * r.co2);
    eqn(4) = (-C(5)) + (tau * r.dmc); 
    eqn(5) = (-C(2)) + (tau * r.eg);
    eqn(6) = (-C(6)) + (tau * r.me);
end
```

## The Recycle Stream 

To recycle stream flowrate is one of the last things to calculate, know that we know the effluence concentrations and flowrates.

[get_recycle_and_fresh_flowrates](https://github.com/wesleyZero/capstone_II/blob/main/reactor_fxns.m#L85-L111)
```matlab
function [F_fresh, R] = get_recycle_and_fresh_flowrates(F_virt_feed, F_real_effluent)
    flow_fxns = flowrate_fxns();

    % Initialize 
    R = flow_fxns.get_blank_flowstream();
    F_fresh = flow_fxns.get_blank_flowstream();

    % Recycle flow
    R.ethylene_carbonate.mol = F_real_effluent.ethylene_carbonate.mol;
    R.methanol.mol = F_real_effluent.methanol.mol;
    R.carbon_dioxide.mol = F_real_effluent.carbon_dioxide.mol;
    R = flow_fxns.set_F_mol(R);

    % Fresh feed flow 
    F_fresh.ethylene_oxide.mol = F_virt_feed.ethylene_oxide.mol - R.ethylene_carbonate.mol;
    F_fresh.methanol.mol = F_virt_feed.methanol.mol - R.methanol.mol;
    F_fresh.carbon_dioxide.mol = F_virt_feed.carbon_dioxide.mol - R.carbon_dioxide.mol;
    F_fresh = flow_fxns.set_F_mol(F_fresh);

end
```


# Genenate Plots to Answer the Original Questions

**functions here are linked, but not shown directly because they are soo long** 

we use all of those functions above, calling them many times for the [isothermal](https://github.com/wesleyZero/capstone_II/blob/main/run_dmc.m#L90-L196) case and the [isobaric](https://github.com/wesleyZero/capstone_II/blob/main/run_dmc.m#L212-L264) cases. Then we evaluate the economic value of these processes(at all tau, T, P) with the [NPV function](https://github.com/wesleyZero/capstone_II/blob/main/economic_fxns.m#L492-L686) to generate the plots we want, and answer the quesions that we want to answer. 

<table>
  <tr>
    <td align="center"><img src="https://github.com/wesleyZero/capstone_II/blob/main/readme/img/ReactorInlet.png" ><br><em>Flowrates going into the plant</em></td>
    <td align="center"><img src="https://github.com/wesleyZero/capstone_II/blob/main/readme/img/SConversion.png" ><br><em>How selective is our primary reaction compared to side reactions</em></td>
  </tr>
  <tr>
    <td align="center"><img src="https://github.com/wesleyZero/capstone_II/blob/main/readme/img/isothermal_F_effluent_74Bar.png" ><br><em>Flowrates coming out of the reactor</em></td>
    <td align="center"><img src="https://github.com/wesleyZero/capstone_II/blob/main/readme/img/isothermal_F_rxtr_total74Bar.png" ><br><em>Total flowrate going into the reactor</em></td>
  </tr>
  <tr>
    <td align="center"><img src="https://github.com/wesleyZero/capstone_II/blob/main/readme/img/isothermal_NPV_all_pressures.png" ><br><em>Econonmic value of different pressures</em></td>
    <td align="center"><img src="https://github.com/wesleyZero/capstone_II/blob/main/readme/img/isothermal_V_reactor.png"><br><em>Reactor volumes at different pressures</em></td>
  </tr>
  <tr>
    <td align="center"><img src="https://github.com/wesleyZero/capstone_II/blob/main/readme/img/separation_feed_composition.png" ><br><em>Composition of whats going into our separation system</em></td>
    <td align="center"><img src="https://github.com/wesleyZero/capstone_II/blob/main/readme/img/isothermal_recycle_74Bar.png"><br><em>Flowrates of our recycle stream</em></td>
  </tr>
</table>

# Nomenclature

## Supercritical_fluids
<img src="https://github.com/wesleyZero/capstone_II/blob/main/readme/img/supercritical_phase_diagram.png" width="500">

**A supercritical fluid** is a fluid that behaves like a gas and a liquid simultaneously. Depending on the thermodynamic conditions (the temperature and pressure) it will behave more like one or the other. The details of this are not what I want the focus to be, thats a topic of physics. Whats important is that we have a reactive component that changes density, **effectively changing the concentration of all of the species in the reactor depending on what thermodynamic conditions we specifiy**. For CO2 we get this type of bizarro phase behavior when we get to pressures of about 76 times atmospheric. 

The function to model this is below. The functions you see are empirical regressions of experimental data since there generally isn't good physical models of these exotic phases of matter. 

[get_supercritical_c02_density function](https://github.com/wesleyZero/capstone_II/blob/main/reactor_fxns.m#L368)
```matlab
function rho = get_supercritical_c02_density(T, P, opt)
    % Input: condition = T or P. Depending on option
    %   P [=] bar
    %   T [=] celcius   
    % Assumptions:
    %   Isobaric model is at 150 bar
    %   Isothermal model is at 140 C
    % Ranges of input
    %   P = [50 bar, 150 bar]
    %   T = [80 C, 140 C]
    % Output: 
    %   rho [=] kg / m^3
    withinTempRange = @(T) T >= 80 && T <= 140;
    withinPressureRange = @(P) P >= 50 && P <= 150;
    
    rho.units = 'kg / m^3';
    switch opt
        case 'isothermal'
            if withinPressureRange(P)
                rho = 1.6746 * P - 12.592;
                    % NIST / Excel Regression
            else
                rho = NaN;
                disp("get_supercritical_c02_density : ERROR : P out of range")
            end
        case 'isobaric'
            if withinTempRange(T)
                if P < 125 % [ Bar ]
                    rho = 356.08 * exp(-0.006 * T);
                        % NIST Data at 100 bar
                else
                    rho = 838.87 * exp(-0.009 * T);
                        % NIST Data at 150 bar
                end
            else
                rho = NaN;
                disp("get_supercritical_c02_density : ERROR : T out of range")
            end
        otherwise
            disp("SUPERCRITICAL C02 DENSITY FUNCTION ERROR: invalid opt")
            rho = NaN;
    end
end
```

## Links
[PDF of our final Report](https://github.com/wesleyZero/capstone_II/blob/main/readme/pdf/capstone2_DMC.pdf)
