from bs4 import BeautifulSoup
import sys
fname = sys.argv[1]
html_doc = open(fname).read()
soup = BeautifulSoup(html_doc)

# check if already have google analytics script
result = [i for i in soup.find_all("script") if str(i).find('UA-21871925-1') >= 0]
if len(result) > 0:
    print >> sys.stderr, "already have Google analytics script, skipping..."
    print(str(soup))
    sys.exit(0)

scriptContent = """
   var _gaq = _gaq || [];
  _gaq.push(['_setAccount', 'UA-21871925-1']);
  _gaq.push(['_trackPageview']);

  (function() {
    var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;
    ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';
    var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(ga, s);
  })();
"""
new_tag = soup.new_tag("script")
new_tag["type"] = "text/javascript"
new_tag.string = scriptContent

tag = soup.body
tag.append(new_tag)
print(str(soup))