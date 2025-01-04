"""
moduletest.jl is a part of the Web-Site engine PkXplore

Copyright or © or Copr. UCBL Lyon, France;  
contributor : [Jean-Pierre Flandrois] ([2024/12/20])
[JP.flandrois@univ-lyon1.fr]

This software is a computer program whose purpose is to create the Web-Site PkXplore.

This software is governed by the [CeCILL|CeCILL-B|CeCILL-C] license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the [CeCILL|CeCILL-B|CeCILL-C]
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the [CeCILL|CeCILL-B|CeCILL-C] license and that you accept its terms.

"""

module App
using GenieFramework
using Stipple

@genietools

@app begin
    @in choix = false
    @in Range_Markers_r = RangeData(0:100)
    @out selectionmsainf=0
    @out selectionmsasup=100
    # @out Range_Markers_r.min=0
    # @out Range_Markers_r.max=30
    # @out Range_Markers_labels =
    #     [Dict(:value => i, :label => string(i) * "%") for i = 0:5:30]
    @onbutton choix begin
        selectionmsainf = Range_Markers_r.range.start
        selectionmsasup = Range_Markers_r.range.stop
        println("cool de ",selectionmsainf," a ",selectionmsasup)
    end 
end
# @app begin
#     @in range_data::R{RangeData{Int}} = RangeData(1:100)

#     # @out Range_Markers_labels =
#     #      [Dict(:value => i, :label => string(i) * "%") for i = 0:10:100]
#     # @onchangeany selectionintervalle
#     #     println(selectionintervalle, typeof(selectionintervalle))
#     #  end
#     println(range_data.max)
# end



function ui()
    [
    #Stipple.range(1:1:30, :Range_Markers_r, markers = true, label = true),
    Stipple.range(1:10:100, :Range_Markers_r, var"marker-labels" = true, color = "orange"),
    # Stipple.range(
    #     0:5:30,
    #     :Range_Markers_r,
    #     markers = true,
    #     #var":marker-labels" = "Range_Markers_labels",
    #     color = "secondary",
    # ),
    cell([" --- info ml {{Range_Markers_r}} *** {{selectionmsasup}} {{selectionmsainf}}",btn("Select", @click(:choix))]) #{ "min": 8, "max": 21 }
]
#     # Html.div(cell([h5("Select a zone and redo"), Stipple.range(
#         #     0:5:100,
#         #     :selectionintervalle,
#         #     markers = true,
#         #     var:"marker-labels" = true,
#         #     color = "secondary"
#         #     ) ]), @showif("trim_fait"))
#     [
#         cell([
#         h4("portion")

#         range(1:10:100,
#               :range_data;
#               label=true,
#               labelalways=true,
#             #   labelvalueleft=Symbol("'Min pos: ' + range_data.min"),
#             #   labelvalueright=Symbol("'Max pos: ' + range_data.max"),
#               color = "secondary")
#     ])
#         cell([" --- info ml {{range_data.min}} {{range_data.max}}"]) #"Process si: {{selectionintervalle}}
    
#     ]
    
#     # [row([
#     #     column(class="q-pa-sm",select(:banqueselectionnée, options=:listebanques, emitvalue=true, clearable=true, useinput=true, counter = true, fillinput=true, filled = true, label="Bank Selection"), sm=6),
#     #     column(class="q-pa-sm",cell([p("Nb seqs"),slider(30:5:75, :requestedseq, var"marker-labels" = true, color = "grey")]), sm=5),
#     #     #column(class="q-pa-sm",btn("Select", @click(:choixbanque)), sm=1),
#     #     ])
#     #cell([" --- info ml {{marker-labels}}"]) #"Process si: {{selectionintervalle}}
#     # row([
#     #     column(span("Hello"), sm=11, class = "bg-blue-2"),
#     #     column(span("Genie"), sm=1, class = "bg-red-2"),
#     # ]),
#     # row([
#     #     column(span("Hello"), sm=11, class = "bg-blue-2"),
#     #     column(span("Genie"), sm=1, class = "bg-red-2"),
#     # ]),
#     # item([
#     # itemsection(avatar = "", icon("volume_up", color = "teal")),
#     # itemsection(slider(0:1:10, :SliderIcon_volume, label = "", color = "teal")),
# # ])

end
@page("/moduletest",ui)
end