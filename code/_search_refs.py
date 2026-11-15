import json, sys, urllib.request, urllib.parse
Q = [
    "Villareal Weight loss, exercise, or both and physical function in obese older adults New England Journal of Medicine 2011",
    "Heymsfield Weight loss composition is one-fourth fat-free mass: a critical review and critique of this widely cited rule",
    "Topp Promislow Devries Ludvik model of beta-cell mass, insulin, and glucose kinetics: pathways to diabetes Journal of Theoretical Biology",
    "Kennedy Berger Geroscience: linking aging to chronic disease Cell 2014",
    "Kennedy Berger Miller Aging Cell 2014 Geroscience: linking aging to chronic disease",
    "Kalyani Corriere Ferrucci Age-related and disease-related muscle loss: the effect of diabetes, obesity, and other diseases Lancet Diabetes Endocrinology",
    "Neeland Linge Birkenfeld lean body mass GLP-1 mitigation strategies",
    "Linge Birkenfeld Neeland Muscle mass and glucagon-like peptide-1 receptor agonists: adaptive or maladaptive response to weight loss",
]
import time
for q in Q[1:]:
    time.sleep(4)
    url = "https://api.crossref.org/works?rows=3&select=DOI,title,author,issued,container-title&query.bibliographic=" + urllib.parse.quote(q)
    with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "reference-check/1.0"}), timeout=30) as r:
        items = json.load(r)["message"]["items"]
    print("Q:", q[:80])
    for it in items:
        a = (it.get("author") or [{}])[0].get("family", "")
        print("   ", it["DOI"], "|", a, it["issued"]["date-parts"][0][0], "|", (it.get("title") or [""])[0][:90], "|", (it.get("container-title") or [""])[0][:30])
