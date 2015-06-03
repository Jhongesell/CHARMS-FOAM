# CHARMS-FOAM
=======================================================================
Computational Hydraulics And Research Modeling Suite based on OpenFOAM

-----------------------------------------------------------------------
0. This library converts the raw markup to HTML. See the list of [supported markup formats](#markups) below.
0. The HTML is sanitized, aggressively removing things that could harm you and your kin—such as `script` tags, inline-styles, and `class` or `id` attributes. See the [sanitization filter](https://github.com/jch/html-pipeline/blob/master/lib/html/pipeline/sanitization_filter.rb) for the full whitelist.
0. Syntax highlighting is performed on code blocks. See [github/linguist](https://github.com/github/linguist#syntax-highlighting) for more information about syntax highlighting.
0. The HTML is passed through other filters in the [html-pipeline](https://github.com/jch/html-pipeline) that add special sauce, such as [emoji](https://github.com/jch/html-pipeline/blob/master/lib/html/pipeline/emoji_filter.rb), [task lists](https://github.com/github/task_list/blob/master/lib/task_list/filter.rb), [named anchors](https://github.com/jch/html-pipeline/blob/master/lib/html/pipeline/toc_filter.rb), [CDN caching for images](https://github.com/jch/html-pipeline/blob/master/lib/html/pipeline/camo_filter.rb), and  [autolinking](https://github.com/jch/html-pipeline/blob/master/lib/html/pipeline/autolink_filter.rb).
0. The resulting HTML is rendered on GitHub.com.
