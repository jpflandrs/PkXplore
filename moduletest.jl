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