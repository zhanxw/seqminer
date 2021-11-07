from bs4 import BeautifulSoup
import sys
fname = sys.argv[1]
html_doc = open(fname).read()
soup = BeautifulSoup(html_doc)

# check if already have google analytics script
result = [i for i in soup.find_all("script") if str(i).find('toc.min.js') >= 0]
if len(result) > 0:
    print >> sys.stderr, "already have toc.min.js script, skipping..."
    print(str(soup))
    sys.exit(0)

result = [i for i in soup.find_all("script") if str(i).find('jquery.min.js') >= 0]
hasJquery = False
if len(result) > 0:
    print >> sys.stderr, "has jquery"
    hasJquery = True

## add css
new_tag = soup.new_tag("link")
new_tag["rel"] = "stylesheet"
new_tag["href"] = "css/toc.css"
tag = soup.head
tag.append(new_tag)

## import jquery
if not hasJquery:
  new_tag = soup.new_tag("script")
  new_tag["src"] = "http://ajax.googleapis.com/ajax/libs/jquery/2.1.1/jquery.min.js"
  tag = soup.head
  tag.append(new_tag)

## import toc.js
new_tag = soup.new_tag("script")
new_tag["src"] = "js/toc.min.js"
tag = soup.head
tag.append(new_tag)

## add <div id = "toc">
new_tag = soup.new_tag("div")
new_tag["id"] = "toc"
tag = soup.body
tag.append(new_tag)

## enable toc.js
new_tag = soup.new_tag("script")
new_tag["type"] = "text/javascript"
new_tag.string = "$('#toc').toc();"
tag = soup.body
tag.append(new_tag)

print(str(soup))