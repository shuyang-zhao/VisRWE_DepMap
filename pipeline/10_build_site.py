"""Step 10 - assemble the self-contained explorer page for GitHub Pages: docs/index.html.

Inlines web/style.css, web/i18n.js, web/app.js, the vendored vis-network build and the data bundle of step 9 into one HTML file
(works offline and from file://).
"""
import os

from _common import banner, paths

banner('10. build docs/index.html')
read = lambda *p: open(os.path.join(*p), encoding='utf-8').read()  # noqa: E731
page = read(paths.WEB, 'index.html')
parts = {
    '%%STYLE%%': read(paths.WEB, 'style.css'),
    '%%I18N%%': read(paths.WEB, 'i18n.js'),
    '%%APP%%': read(paths.WEB, 'app.js'),
    '%%DATA%%': read(paths.INTERMEDIATE, 'explorer_data.json'),
    '%%VIS%%': read(paths.WEB, 'vendor', 'vis-network-9.1.9.min.js'),
}
for k, v in parts.items():                      # vendor code last, so no placeholder can appear inside replaced content
    assert page.count(k) == 1, k
for k in ('%%STYLE%%', '%%I18N%%', '%%APP%%', '%%DATA%%', '%%VIS%%'):
    page = page.replace(k, parts[k].replace('</script>', '<\\/script>') if k != '%%STYLE%%' else parts[k])
out = os.path.join(paths.ROOT, 'docs', 'index.html')
os.makedirs(os.path.dirname(out), exist_ok=True)
open(out, 'w', encoding='utf-8').write(page)
open(os.path.join(paths.ROOT, 'docs', '.nojekyll'), 'w').close()
print(f'wrote docs/index.html ({os.path.getsize(out) / 1e6:.1f} MB)')
