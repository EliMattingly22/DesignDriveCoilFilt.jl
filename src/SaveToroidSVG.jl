function Coords2StringList(R,Z;SigDigits = 4)
    L = length(R)
    Str = " "
    for i in 1:L
        Str = Str*"$(round(R[i];sigdigits=SigDigits)), $(round(Z[i];sigdigits=SigDigits)) "
    end
    return Str
end

function DCoreParams2SVG(ID,OD,FileName;NPts = 100)
    r, z =DCoreGeom(ID/2, OD/2;NPts = NPts,PlotOn=false,UpsamplePoints = 1e4)
    z = z .+ maximum(z)
    Str = Coords2StringList(r,z;SigDigits = 4)
    open(FileName,"w") do file
        write(file,
"""
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg
   xmlns:dc="http://purl.org/dc/elements/1.1/"
   xmlns:cc="http://creativecommons.org/ns#"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
   xmlns:svg="http://www.w3.org/2000/svg"
   xmlns="http://www.w3.org/2000/svg"
   xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"
   xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
   width="$(round(OD/2))mm"
   height="$(round(maximum(abs.(z))))mm"
   viewBox="0 0 $(round(OD/2)) $(round(maximum(abs.(z))))"
   version="1.1"
   id="svg1501"
   inkscape:version="1.0.2 (e86c870879, 2021-01-15, custom)"
   sodipodi:docname="drawing.svg">
  <defs
     id="defs1495" />
  <sodipodi:namedview
     id="base"
     pagecolor="#ffffff"
     bordercolor="#666666"
     borderopacity="1.0"
     inkscape:pageopacity="0.0"
     inkscape:pageshadow="2"
     inkscape:zoom="0.35"
     inkscape:cx="400"
     inkscape:cy="560"
     inkscape:document-units="mm"
     inkscape:current-layer="layer1"
     inkscape:document-rotation="0"
     showgrid="false"
     inkscape:window-width="1624"
     inkscape:window-height="1228"
     inkscape:window-x="190"
     inkscape:window-y="190"
     inkscape:window-maximized="0" />
  <metadata
     id="metadata1498">
    <rdf:RDF>
      <cc:Work
         rdf:about="">
        <dc:format>image/svg+xml</dc:format>
        <dc:type
           rdf:resource="http://purl.org/dc/dcmitype/StillImage" />
        <dc:title></dc:title>
      </cc:Work>
    </rdf:RDF>
  </metadata>
  <g
     inkscape:label="Layer 1"
     inkscape:groupmode="layer"
     id="layer1">
    <path
       style="fill:none;stroke:#000000;stroke-width:0.5px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;image-rendering:auto"
       d="M  $(Str) z"
       id="path2064" />
  </g>
</svg> """)
    end
end


