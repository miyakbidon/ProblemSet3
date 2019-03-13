# #Definitions of metabolites
# "1. Aspartate"
# "2. Citrulline"
# "3. 2-(Nomega-arginino)succinate"
# "4. Fumarate"
# "5. Arginine"
# "6. Ornithine"
# "7. Urea"
# "8. Carbamoyl phosphate"
# "9. O2"
# "10. NADPH"
# "11. NADP+"
# "12. H+"
# "13. H2O"
# "14. NO"
# "15. ATP"
# "16. AMP"
# "17. PPi"
# "18. Pi"
#
# #Reactions
# "V1: ATP + Citrulline + Aspartate -> AMP + PPi + 2-(Nomega-L-arginine)succinate"
# "V2: 2-(Nomega-L-arginine)succinate -> Fumarate + Arginine"
# "V3: Arginine + H2O -> Ornithine + Urea"
# "V4: Carbamoyl phosphate + Ornithine -> Orthophosphate + Citrulline"
# "V5: 2 Arginine + 4 O2 + 3 NADPH + 3 H+ <-> 2 NO + 2 Citruilline + 3 NADP+ + 4 H2O"
# "B1: Carbamoyl phosphate (b) -> Carbamoyl phosphate"
# "B2: Aspartate (b) -> Aspartate"
# "B3: Fumarate -> Fumarate (b)"
# "B4: Urea -> Urea (b)"
# "B5: ATP (b) -> ATP"
# "B6: AMP -> AMP (b)"
# "B7: PPi -> PPi (b)"
# "B8: Pi -> Pi (b)"
# "B9: O2 (b) -> O2"
# "B10: NADPH (b) -> NADPH"
# "B11: H+ (b) -> H+"
# "B12: NO -> NO (b)"
# "B13: H2O -> H2O (b)"       # from V5
# "B14: NADP+ -> NADP+ (b)"
# "B15: H2O (b) -> H2O"       # from V3

include("Flux.jl")



# Stoichiometric Matrix
S = Array{Float64}([
    -1.0     0     0     0      0     0    1.0    0     0     0     0     0     0     0     0     0     0     0     0     0
    -1.0     0     0    1.0    2.0    0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     1.0   -1.0    0     0      0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
      0     1.0    0     0      0     0     0   -1.0    0     0     0     0     0     0     0     0     0     0     0     0
      0     1.0  -1.0    0    -2.0    0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
      0      0    1.0  -1.0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
      0      0    1.0    0      0     0     0     0   -1.0    0     0     0     0     0     0     0     0     0     0     0
      0      0     0   -1.0     0    1.0    0     0     0     0     0     0     0     0     0     0     0     0     0     0
      0      0     0     0    -4.0    0     0     0     0     0     0     0     0    1.0    0     0     0     0     0     0
      0      0     0     0    -3.0    0     0     0     0     0     0     0     0     0    1.0    0     0     0     0     0
      0      0     0     0     3.0    0     0     0     0     0     0     0     0     0     0     0     0     0   -1.0    0
      0      0     0     0    -3.0    0     0     0     0     0     0     0     0     0     0    1.0    0     0     0     0
      0      0   -1.0    0     4.0    0     0     0     0     0     0     0     0     0     0     0     0   -1.0    0    1.0
      0      0     0     0     2.0    0     0     0     0     0     0     0     0     0     0     0   -1.0    0     0     0
    -1.0     0     0     0      0     0     0     0     0    1.0    0     0     0     0     0     0     0     0     0     0
     1.0     0     0     0      0     0     0     0     0     0   -1.0    0     0     0     0     0     0     0     0     0
     1.0     0     0     0      0     0     0     0     0     0     0   -1.0    0     0     0     0     0     0     0     0
      0      0     0    1.0     0     0     0     0     0     0     0     0   -1.0    0     0     0     0     0     0     0

]);

# Atom Matrix
# This matrix is arranged in columns as stated in the metabolite definition and in rows as C H N O P S.
A = [
    4    6   10    4    6    5    1    1    0   21   21    0    0    0   10   10    0    0
    7   13   18    4   14   12    4    4    0   30   29    1    2    0   16   14    4    3
    1    3    4    0    4    2    2    1    0    7    7    0    0    1    5    5    0    0
    4    3    6    4    2    2    1    5    2   17   17    0    1    1   13    7    7    4
    0    0    0    0    0    0    0    1    0    3    3    0    0    0    3    1    2    1
    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
];

Check = A * S

# Species Balance
species_bounds_array = [
         0.0    0.0	;       # 1.
         0.0    0.0	;       # 2.
         0.0    0.0	;       # 3.
         0.0    0.0	;       # 4.
         0.0    0.0	;       # 5.
         0.0    0.0	;       # 6.
         0.0    0.0	;       # 7.
         0.0    0.0	;       # 8.
         0.0    0.0	;       # 9.
         0.0    0.0	;       # 10.
         0.0    0.0	;       # 11.
         0.0    0.0	;       # 12.
         0.0    0.0	;       # 13.
         0.0    0.0	;       # 14.
         0.0    0.0	;       # 15.
         0.0    0.0	;       # 16.
         0.0    0.0	;       # 17.
         0.0    0.0	;       # 18.
    ];

# Objective coefficient array
objective_coefficient_array = [
    	0.0    ;    # V1
        0.0    ;    # V2
        0.0    ;    # V3
        0.0    ;    # V4
        0.0    ;    # V5
        0.0    ;    # B1
        0.0    ;    # B2
        0.0    ;    # B3
        1.0    ;    # B4: maximizing Urea production
        0.0    ;    # B5
        0.0    ;    # B6
        0.0    ;    # B7
        0.0    ;    # B8
        0.0    ;    # B9
        0.0    ;    # B10
        0.0    ;    # B11
        0.0    ;    # B12
        0.0    ;    # B13
        0.0    ;    # B14
        0.0    ;    # B15
    ];

# Defining V_max variables for default bounds array
V_max1 = 7.308    ;  # mmol/g-DW*hr
V_max2 = 1.242    ;  # mmol/g-DW*hr
V_max3 = 8.964    ;  # mmol/g-DW*hr
V_max4 = 3.1716   ;  # mmol/g-DW*hr
V_max5 = 0.49319  ;  # mmol/g-DW*hr

# Defining Control variables for default bounds array for the 5 reactions
Cntrl1 = 0.9225 * 0.9898  ;
Cntrl2 = 1                ;
Cntrl3 = 0.1418           ;
Cntrl4 = 0.7373           ;
Cntrl5 = 0.9865           ;

# Default bounds array
default_bounds_array = [
       0	V_max1*Cntrl1	;	# V1
       0	V_max2*Cntrl2	;	# V2
       0	V_max3*Cntrl3	;	# V3
       0	V_max4*Cntrl4	;	# V4
       -V_max5*Cntrl5	V_max5*Cntrl5	;	# V5
       0	10     ;       # B1
       0	10     ;       # B2
       0	10     ;       # B3
       0	10     ;       # B4
       0	10     ;       # B5
       0	10     ;       # B6
       0	10     ;       # B7
       0	10     ;       # B8
     -10	10     ;       # B9
     -10	10     ;       # B10
     -10	10     ;       # B11
     -10	10     ;       # B12
     -10	10     ;       # B13
     -10	10     ;       # B14
       0	10     ;       # B15
   ];

(objective_value, calculated_flux_array, dual_value_array, uptake_array, exit_flag, status_flag) = calculte_optimal_flux_distribution(S,default_bounds_array,species_bounds_array,objective_coefficient_array,min_flag=false)
