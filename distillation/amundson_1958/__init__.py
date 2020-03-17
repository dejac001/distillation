"""
.. graphviz::

    digraph {
        input -> initial_guess -> calc_K -> mass_balance -> bubble_pt -> converg_1;
        converg_1 -> calc_K[label="No",color=red];
        converg_1 -> energy_balance[label="Yes",color=green];
        energy_balance -> converg_2;
        converg_2 -> finish[label="Yes",color=green];
        converg_2 -> mass_balance[label="No",color=red];

        input [shape=polygon,side=4,label="Input equilibrium and enthalpy data.\n Input design variables"];
        initial_guess [shape=polygon,side=4,label="Pick L, V, T\n on every stage"];
        calc_K [shape=polygon,side=4,label="Calculate K for each component\n on each stage"];
        mass_balance [shape=polygon,side=4,label="Solve component mass\n balances in matrix form",color=lightblue,style=filled];
        bubble_pt [shape=polygon,side=4,label="Calculate T on each stage\n using bubble point calculation"];
        converg_1 [shape=diamond,label="Determine whether\n T is converged\n for all stages"];
        energy_balance [shape=polygon,side=4,label="Use energy balance to calculate\n L, V on every stage",color=lightblue,style=filled];
        converg_2 [shape=diamond,label="Determine whether\n L, V are converged\n for all stages"];
        finish [shape=ellipse,label="Finished"];
    }
"""
