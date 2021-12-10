'''
Paginator
Author: Jake Bradford

- Parses an iterable as pages
- Useful for partitioning iterators containing many items
- Can be used on any object that it is iterable
- Page length and start page can be specified

With thanks to RobertB (8 September 2017) for the inspiration:
    https://stackoverflow.com/a/46107096
'''

#import random
#import string
#random.seed(20210209)

class Paginator:
    def __init__(self, iterable, page_len, start_page=0):
        self.iterable = iterable
        self.page_len = page_len
        self.current_page = 0
        self.desired_page = start_page
        self.current_item = 0

    def __iter__(self):
        page = []
        
        if self.page_len == 0:
            yield 1, self.iterable
        else:
            for i in self.iterable:
                self.current_item += 1
                
                if self.current_page < self.desired_page:
                    if self.current_item % self.page_len == 0:
                        self.current_page += 1
                    continue
            
                page.append(i)
                if len(page) == self.page_len:
                    yield self.current_page, page 
                    page = []
                    self.current_page += 1
            if page:
                yield self.current_page, page
                self.current_page += 1


#a = {}
#for b in random.sample(range(1, 1000000), 30):
#    a[b] = b*2
#
#
#for idx, i in Paginator(a, 7, 2):
#    print(idx)
#    for j in i:
#        print(f'\t{j}')
