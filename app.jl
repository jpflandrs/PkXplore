"""
app.jl is a part of the Web-Site engine PkXplore

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
using Main.Analysis

using GenieFramework.Genie.Requests: postpayload

route("/") do
    [ h4("Vous voici sur le serveur PkXplorer") ]
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


include("nucworkshop.jl")
include("moduletest.jl")
end

