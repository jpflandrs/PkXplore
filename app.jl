module App
using GenieFramework
using Main.Analysis

using GenieFramework.Genie.Requests: postpayload

route("/") do
    [ h4("Welcome to the fasta analysis module"),
      a("fasta form", href="/form") ]
end

route("/form") do
    Html.form(action = "/result", method="POST", [
        input(type="textarea",rows="10", cols="150", name="N", placeholder="Fasta")
        input(type="submit", value="Send")
    ])
#     html("""
#   <h1>FASTA Sequence Processor</h1>
#   <form action="/blast" method="POST">
#       <label for="fasta_text">Enter FASTA sequence:</label><br>
#       <textarea name="fasta_text" rows="10" cols="50"></textarea><br>
#       <br>
#       <input type="submit" value="Process">
#   </form>
#   """)
end


route("/result", method=POST) do
    N = postpayload(:N)
    
    p("$N", style="font-size:20px")
end


include("ReactiveForm.jl")
include("moduletest.jl")
end

