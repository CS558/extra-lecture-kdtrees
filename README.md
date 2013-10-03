CS558: kd trees
===============

# Problem statement

An orthogonal range, $I \subset \mathbb{R}^d$, is a product of (possibly unbounded) intervals.  Let $S$ be a set of $n$ points in Euclidean space.  Then the orthogonal range searching problem is to find the set of all points in $S$ which are contained in $I$.  There are a few variations on how we could process the result here.  We might count all the points, compute some aggregate quantity (ie map/reduce type of operation), or even just reporting all the points that happen to lie in the range.  We will focus on the latter version

There is an obvious O(n) algorithm to solve this problem, and if we don't do anything to preprocess the points this is obviously optimal.  In fact, because any search could possibly return all the points in an interval, it would seem that O(n) is actually a tight bound on the worst case performance.  However, if we are a bit more careful we can separate out the complexity as a function of both the number of points in the data structure, as well as the number of points that are returned (which we shall call k).  Under this more detailed analysis, the complexity of naive range searching is O(n + k).

The goal of this lecture and the next one is to look at some ways to speed up orthogonal range queries using precomputation.


## 1D range searching

To get started, consider the special case where d=1.  In this case, there is an obvious solution:  just sort the data and use binary search.  This takes O(n) space, requires O(n log(n)) preprocessing time, and answers queries in O(log(n) + k) time.

## Higher dimensions

It would be nice to extend this to higher dimensions, but it is a fact that there is no complete total order on R^d for d>1.  As a result, binary search doesn't directly generalize to higher dimensions, and so we need to use some other data structure.  It is also not possible in general to get O(log(n)) query overhead with O(n) storage in higher dimensions.  The best results known are as follows:

* O(n) storage, O(n^{1-1/d} + k) query
* O(n log^d(n)) storage, O(log^d n) query

This is basically optimal (within a log factor on the latter version), as proven by Chazelle.  Both of these results are also acheivable in practice.  Today we are going to talk about the first case, which is a bit simpler and can be achieved using kd-trees

# kd-trees

The k in kd-trees originally stood for the of them, and they were invented by Bentley in the 1970s as a way to generalize binary search trees.  A kd-tree is a balanced binary tree, where each node in the tree stores a point along with pointers to the left and right children.  Each node in a kd-tree must satisfy the following properties:

* The tree is balanced, ie number of nodes in left subtree = number of nodes in right subtree +/- 1
* If and only if a point's d'th component is less than p's then it is in the left subtree.
* d(n) = (d(p(n)) + 1) % dimension

Draw a picture

## Construction
Using the above invariant, it is pretty simple to see how we could construct such a data structure.  Just recursively split the point set by the median, divide into two equal sized subsets, and repeat:

```javascript
function makeTree(d, points) {
  if(points.length === 0) {
    return null
  }
  points.sort(function(a,b) {
    return a[d] - b[d]
  })
  var mid = points.length >> 1
  return {
    d: d,
    p: points[mid],
    left: makeTree((d+1)%dimension, points.slice(0, mid)),
    right: makeTree((d+1)%dimension, points.slice(mid+1))
  }
}
```

The above solution works, and you can easily verify that it satisifes the conditions we placed on the kd-tree.  But we can speed it up a bit, and the trick here is how we select the midpoint.  Using the naive sort-and-split method, we do O(n log(n)) work per level, giving a total complexity of O(n log^2 (n)).  We can do a bit better though if we presort all the points along each axis and maintain d sorted lists, giving an O(d n log(n)) algorithm.  However, the best solution is really to use a fast median selection algorithm to get the split point, which takes O(n) work per level giving an overall O(n log(n)) algorithm.

## Range query

The key to performing a range query on a kdtree is to understand how a kd-tree partitions space.  The basic idea is that we take our range, insert it into the tree, split it along the plane into two nonoverlapping subranges and recurse.  If we ever get to a node where the range does not intersect the bounds of that particular subtree, then we can stop early and avoid wasting time visiting extra nodes.  In JavaScript, this is what it would look like:

```javascript
function queryBox(tree, lo, hi) {
  //Check if terminal node
  if(!tree) {
    return []
  }
  //Check if range is emtpy
  for(var i=0; i<dimension; ++i) {
    if(lo[i] > hi[i]) {
      return []
    }
  }
  //Check if point is contained in range
  var p = tree.p
  var result = [ p ]
  for(var i=0; i<dimension; ++i) {
    if(lo[i] > p[i] ||  p[i] > hi[i]) {
      result = []
    }
  }
  //Split range along d
  var d = tree.d
  var left_hi = hi.slice()
  var right_lo = lo.slice()
  left_hi[d] = Math.min(p[d], hi[d])
  right_lo[d] = Math.max(p[d], lo[d])
  //Recursively visit subtrees
  return result.concat(queryBox(tree.left, lo, left_hi)).concat(queryBox(tree.right, right_lo, hi))
}
```

## Analysis
The question then is how much overhead this takes.  In the best possible case, the range is completely contained in just one of the subtrees, and so we could at most answer the entire range query in no less than Omega(log(n) + k) time since the tree is balanced.  The more interesting question though is what happens in the worst case.  To see what happens here, consider querying against a very long and thin range.  Here the overhead is going to be proportional to the number of lines that the range intersects.  Taking the limit where the range reduces to a line, let us try to bound the number of intersections of a line with the partitions of a kdtree.

Starting from a root node, let us call the number of crossings C(n).  Then, it is true that:

C(n) = O(1) + 2 * C(n/4)

Which by the master theorem, it is true that:

C(n) = O(sqrt(n))

So the overall cost of a range query is O(sqrt(n) + k) in the worst case.  We could extend this same analysis to higher dimension to see that,

C(n) = O(1) + 2^(d-1) C(n / 2^d)

And so in general, the cost of range searching on a kdtree is:

O(n^(1 - 1/d) + k)

As a result, if d becomes large then the speed up from kd-trees tends to 0 quickly.  This is actually optimal for any data structure that uses at most O(n) extra space.  One could also ask how frequently the bad cases pop up.  On the other hand, if the size of the query regions is relatively small, then we would expect to get closer to O(log(n)) time for most of the queries.

## Advanced stuff

In addition to range queries, there are a few other things that we can do with kd-trees.

### Closest point
Probably the number one application for kd-trees is nearest neighbor queries.  The basic idea is that we maintain and iteratively update a bounding sphere around a point as we traverse the tree, visiting only the nodes which are close to the point.

