# -*- coding: utf-8 -*-
"""Fetch every cited DOI from Crossref and check first author, year and title keywords against the intended reference.
Writes results/references_crossref.json (the bibliographic data used by the manuscript) and prints a report.
A reference is used in the manuscript only if its record is found and every check passes."""
import json
import os
import re
import time
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.abspath(os.path.join(HERE, "..", "results"))

# key: (doi, first-author surname, year, title keywords)
CAND = {
    "Yang2023": ("10.1016/j.cell.2022.12.027", "Yang", 2023, ["epigenetic information", "mammalian aging"]),
    "Levine2018": ("10.18632/aging.101414", "Levine", 2018, ["epigenetic biomarker", "lifespan"]),
    "Horvath2013": ("10.1186/gb-2013-14-10-r115", "Horvath", 2013, ["methylation age", "tissues"]),
    "Hannum2013": ("10.1016/j.molcel.2012.10.016", "Hannum", 2013, ["methylation profiles", "aging"]),
    "Lu2023": ("10.1038/s43587-023-00462-6", "Lu", 2023, ["methylation age", "mammalian"]),
    "Bell2019": ("10.1186/s13059-019-1824-y", "Bell", 2019, ["methylation", "clocks"]),
    "LopezOtin2023": ("10.1016/j.cell.2022.11.001", "López-Otín", 2023, ["hallmarks", "aging"]),
    "Wilding2021": ("10.1056/NEJMoa2032183", "Wilding", 2021, ["semaglutide", "obesity"]),
    "Jastreboff2022": ("10.1056/NEJMoa2206038", "Jastreboff", 2022, ["tirzepatide", "obesity"]),
    "Sargeant2019": ("10.3803/EnM.2019.34.3.247", "Sargeant", 2019, ["lean body mass", "glucagon"]),
    "Neeland2024": ("10.1111/dom.15728", "Neeland", 2024, ["lean body mass", "glucagon"]),
    "Heymsfield2014": ("10.1111/obr.12143", "Heymsfield", 2014, ["weight loss composition", "fat-free"]),
    "Cruz2019": ("10.1093/ageing/afy169", "Cruz-Jentoft", 2019, ["sarcopenia", "consensus"]),
    "Srikanthan2014": ("10.1016/j.amjmed.2014.02.007", "Srikanthan", 2014, ["muscle mass index", "longevity"]),
    "Srikanthan2010": ("10.1371/journal.pone.0010805", "Srikanthan", 2010, ["sarcopenia", "insulin resistance"]),
    "DeFronzo2009": ("10.2337/dc09-S302", "DeFronzo", 2009, ["skeletal muscle insulin resistance", "primary defect"]),
    "Linge2024": ("10.1161/circulationaha.124.067676", "Linge", 2024, ["muscle mass", "glucagon"]),
    "Wolfe2006": ("10.1093/ajcn/84.3.475", "Wolfe", 2006, ["muscle", "health and disease"]),
    "Brownlee2001": ("10.1038/414813a", "Brownlee", 2001, ["biochemistry", "diabetic complications"]),
    "Hotamisligil2006": ("10.1038/nature05485", "Hotamisligil", 2006, ["inflammation", "metabolic disorders"]),
    "Tchkonia2010": ("10.1111/j.1474-9726.2010.00608.x", "Tchkonia", 2010, ["fat tissue", "senescence"]),
    "Batsis2018": ("10.1038/s41574-018-0062-9", "Batsis", 2018, ["sarcopenic obesity", "older adults"]),
    "Villareal2011": ("10.1056/nejmoa1008234", "Villareal", 2011, ["weight loss", "exercise", "obese older adults"]),
    "Heymsfield2021": ("10.1001/jamanetworkopen.2020.33457", "Heymsfield", 2021, ["bimagrumab", "fat mass"]),
    "Barzilai2016": ("10.1016/j.cmet.2016.05.011", "Barzilai", 2016, ["metformin", "aging"]),
    "Simental2008": ("10.1089/met.2008.0034", "Simental", 2008, ["triglycerides", "glucose", "insulin resistance"]),
    "Zhao2023": ("10.1038/s41598-023-43512-0", "Zhao", 2023, ["triglyceride", "glucose"]),
    "Bergman1979": ("10.1152/ajpendo.1979.236.6.E667", "Bergman", 1979, ["insulin sensitivity"]),
    "Topp2000": ("10.1006/jtbi.2000.2150", "Topp", 2000, ["cell mass", "glucose kinetics"]),
    "Sobol2001": ("10.1016/S0378-4754(00)00270-6", "Sobol", 2001, ["global sensitivity indices"]),
    "Saltelli2010": ("10.1016/j.cpc.2009.09.018", "Saltelli", 2010, ["variance based sensitivity"]),
    "Jansen1999": ("10.1016/S0010-4655(98)00154-4", "Jansen", 1999, ["analysis of variance designs"]),
    "Marino2008": ("10.1016/j.jtbi.2008.04.011", "Marino", 2008, ["uncertainty and sensitivity", "systems biology"]),
    "Raue2009": ("10.1093/bioinformatics/btp358", "Raue", 2009, ["identifiability", "profile likelihood"]),
    "Petzold1983": ("10.1137/0904010", "Petzold", 1983, ["stiff and nonstiff"]),
    "Virtanen2020": ("10.1038/s41592-019-0686-2", "Virtanen", 2020, ["scipy"]),
    "Harris2020": ("10.1038/s41586-020-2649-2", "Harris", 2020, ["array programming", "numpy"]),
    "Hunter2007": ("10.1109/MCSE.2007.55", "Hunter", 2007, ["matplotlib"]),
    "Kirkwood1977": ("10.1038/270301a0", "Kirkwood", 1977, ["evolution of ageing"]),
    "Hall2011": ("10.1016/S0140-6736(11)60812-X", "Hall", 2011, ["energy imbalance", "bodyweight"]),
    "Lu2020": ("10.1038/s41586-020-2975-4", "Lu", 2020, ["reprogramming", "youthful epigenetic"]),
    "Chan2025": ("10.18632/aging.206220", "Chan", 2025, ["methylation entropy"]),
    "Fahy2019": ("10.1111/acel.13028", "Fahy", 2019, ["reversal of epigenetic aging"]),
    "Kennedy2014": ("10.1016/j.cell.2014.10.039", "Kennedy", 2014, ["geroscience", "chronic disease"]),
    "Prado2024": ("10.1016/S2213-8587(24)00272-9", "Prado", 2024, ["muscle matters"]),
    "Kalyani2014": ("10.1016/s2213-8587(14)70034-8", "Kalyani", 2014, ["muscle loss", "diabetes", "obesity"]),
}


def norm(s):
    return re.sub(r"[^a-z0-9 ]", "", s.lower().replace("-", " ").replace("‐", " "))


def main():
    cache_path = os.path.join(RES, "references_crossref.json")
    cache = json.load(open(cache_path, encoding="utf-8")) if os.path.exists(cache_path) else {}
    for key, (doi, sur, year, kws) in CAND.items():
        if key in cache and cache[key].get("status") == "ok":
            continue
        url = "https://api.crossref.org/works/" + urllib.request.quote(doi, safe="/")
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "reference-check/1.0"}), timeout=30) as r:
                msg = json.load(r)["message"]
        except Exception as e:
            cache[key] = dict(status="NOTFOUND", doi=doi, error=str(e)[:80])
            print("%-16s NOT FOUND (%s)" % (key, str(e)[:60]))
            continue
        title = (msg.get("title") or [""])[0]
        first = (msg.get("author") or [{}])[0].get("family", "")
        issued = (msg.get("issued", {}).get("date-parts") or [[None]])[0][0]
        pub = (msg.get("published-print", {}).get("date-parts") or msg.get("published", {}).get("date-parts") or [[None]])[0][0]
        ok_author = norm(sur) in norm(first) or norm(first) in norm(sur)
        ok_year = year in (issued, pub)
        ok_title = all(norm(k) in norm(title) for k in kws)
        status = "ok" if (ok_author and ok_year and ok_title) else "MISMATCH"
        cache[key] = dict(status=status, doi=doi, title=title, authors=[(a.get("family", ""), a.get("given", "")) for a in msg.get("author", [])],
                          journal=(msg.get("container-title") or [""])[0], volume=msg.get("volume"), issue=msg.get("issue"), page=msg.get("page"),
                          article_number=msg.get("article-number"), year=issued, pub_year=pub, type=msg.get("type"))
        print("%-16s %-8s author=%s year=%s/%s | %s" % (key, status, first, issued, pub, title[:80]))
        time.sleep(0.15)
    json.dump(cache, open(cache_path, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    bad = [k for k, v in cache.items() if v["status"] != "ok"]
    print("\nverified: %d, not verified: %s" % (len(cache) - len(bad), bad))


if __name__ == "__main__":
    main()
