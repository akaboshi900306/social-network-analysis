rm(list = ls(all = TRUE))
setwd("C:/Users/dplewi3/Dropbox")
#setwd("~/Dropbox")

# start with base r
system.time(read.csv("C:/Users/akabo/Downloads/reddit_data.csv", header = TRUE))
# almost 90 seconds

# from Hadley Wicham's readr
library(readr)
system.time(read_csv("C:/Users/akabo/Downloads/reddit_data.csv", col_names = TRUE))
# 8 seconds

# data table
library(data.table)
system.time(fread("C:/Users/akabo/Downloads/reddit_data.csv", header = TRUE))
# 3 seconds


# okay, so now let's load this into r and peform a few grouping operations
red = read.csv("C:/Users/akabo/Downloads/reddit_data.csv", header = TRUE)

# first, let's aggregate by user, using base r, to see how many comments they made in each subreddit and when the first time they posted in this subreddit was
system.time(cbind(aggregate(utc ~ username + subreddit, data = red, FUN=length), aggregate(utc ~ username + subreddit, data = red, FUN=min)))
# 65 seconds


# with summaryBy from doBy
library(doBy)
# at least the syntax is a little easier
count_fun = function(x){l=length(x)}
system.time(summaryBy(utc ~ username + subreddit, data = red, FUN=c(min,length)))
# 41 seconds here

# let's try dyplr
library(dplyr)
system.time(red %>%
	group_by(username, subreddit) %>%
	summarise(count = n(), earliest = min(utc)))
# much faster, 2.5 seconds

# data table
red2 = fread("C:/Users/akabo/Downloads/reddit_data.csv", header = TRUE)
system.time(red2[, count := .N, by = c("username", "subreddit")][, min := min(utc), by = c("username", "subreddit")])
# less than a second

# so, fastest on operation, and added bonus:
# with
red2[, count := .N, by = c("username", "subreddit")][, min := min(utc), by = c("username", "subreddit")]
# notice that this has been performed by reference
# i.e., we still have the original data object here intact too, to manipulate as we wish
# if we wanted just the aggregated table, then run the following
# will get into the list notation inside the data table call in the section below
red3 = unique(red2[,list(count = count, min = min), by = c("username", "subreddit")])


# save 100000 random rows of this as a data table to use for later examples
rdt = red3[sample(nrow(red3), 100000, replace = FALSE)]
# and as a data frame
rdf = as.data.frame(rdt)

head(rdf)
rdt

# note also that data table also gets rid of the need for ever calling head()
# automatiicaly gives first and last 5 to save time without cluttering console

# for a data frame, could combine something like this for a workaround
hdtl = function(df){
	rbind(head(df), tail(df))
}

hdtl(rdf)

# can check features of data table
typeof(rdt)
class(rdt)

# data table is a data frame with enhanced features, i.e., still inherits the data frame class

# a lot of the data table commands share some similarity with sql
# so it's helpful to think about the package like base r data frames + some sql-like enhancements

# one main difference is that data frames don't have row names
rownames(rdf) = order(-seq_len(nrow(rdf)))
hdtl(rdf)

# won't throw an error, but will ignore you since this parameter doesn't exist
rownames(rdt) = order(-seq_len(nrow(rdt)))
rdt

# why no row names? 
# use keys instead, keys are like enhanced (silent) row names that can also be applied to groups
setkey(rdt, username)

# automatically sorts on key when set
rdt

# can check a key with 
key(rdt)

# can set multiple keys with
setkeyv(rdt, c("username", "subreddit"))


# keys are incredibly helpful for fast merging operations
# username and subreddit will provide a unique identifying because of our grouping from before

# if we wanted to merge our count and min for these users back into the main data object
setkeyv(red2, c("username", "subreddit"))

# data table knows to automatically merge on the key, without it having to be specified
merge(red2, rdt[,-3], all.x = TRUE)

# quick aside, the use of negative indices here is a "negative match"
# this is performing the duties of our normal indexing in the reverse fashion
# can be helpful for quickly excluding a column without having to manipulate the actual object
# to negative match multiple columns, use parens and sequences, like
red2[,-(2:3)]
# this type of syntax is helpful for excluding the middle columns of an object
# in situations where something like
red2[,1:(ncol(red2)-2)]
# wouldn't work

# back to the merge
# this is actually calling merge.data.tabe since we're providing a dt to the function
# it's worth checking out ?merge.data.table

# all.x in this case is performing a left join, keeping all values of the lefthand object
# even if there's no match in y
# so note we still have 14m rows 

# let's benchmark this
system.time(merge(red2, rdt[,-3], all.x = TRUE))
# less than a second

# compared with the same process using data frames
# here we have to supply the merge parameters direclty
system.time(merge(red, rdf[,-3], by = c("username", "subreddit"), all.x = TRUE))
# about 53 seconds

# the dt version is faster because it can just scan row names
# instead of needing to check each row to see if the condition is true


# you may have noticed some slight differences in our df and dt subsetting operations

# let's say that we have an index we wanted to subset by
# sometimes it's helpful to build a logical index to pass to something else to subset
# useful in cases where we only want specific group members that meet a certain condition
# what if we only wanted users that commented more than once, to get rid of one-time throwaway accounts

index = duplicated(rdt$username)

# for a data frame, we subset using a comma
dim(rdf[index,])

# for a data table, comma is unnecesary
dim(rdt[index])

# let's unpack what this is actually doing, e.g., would it work on columns too
index = rep(c(TRUE, FALSE), 2)
dim(rdf[,index])
dim(rdt[index])
# check the error...

# okay, so the error explains that it's explictly expecting an index on the rows

# also notice the difference here
dim(rdf[-3])
dim(rdt[-3])

# adding some columns to dt and df objects
# for a df, can use the $ operator
# average comments for all users
rdf$avg_comm = mean(rdf$count)
hdtl(rdf)

# for a dt, assignment occurs directly inside of the dt
# when we are adding a new column, we use the := operator
rdt[, avg_comm := mean(count)]
# here, can refer to variable names directly without $

# in a df, we usually refer to columns or rows using statements before and after the , in the [] index
# so to get user subreddit pairs with the most comments, we'd say
rdf[rdf$count == max(rdf$count),]
rdt[count == max(count), most_count := count]
# if doing assignment based on a condition in the rows,
# will get NAs (not 0s) for the rows to which this doesn't apply


# what if we wanted to assign several columns at once?
# just combine bracket operators

# mean, variance, and skewness of counts for all users
# can get skewness from the moments package
library(moments)
rdt[, avg_comm := mean(count)][, var_comm := var(count)][, skew_comm := skewness(count)]

# what if we wanted these numbers to be user-specific? 

# data tables let us add a 3rd argument inside of []
# let's combine some syntax from the grouping operations above
# in the last spot, we add in a "by" indicator as a group marker 
rdt[, avg_comm := mean(count), by = username][, var_comm := var(count), by = username][, skew_comm := skewness(count), by = username]


# what if we wanted to assign a lot of columns at once, according to a function?
# let's take the maximum value of those three summary stats for each subreddit

# first define the column names we will operate on
cols = colnames(rdt)[(ncol(rdt) - 2):ncol(rdt)]

# often the paste function is useful for something like this
max_cols = paste("max_",cols,"_subreddit",sep="")

# so now we have two character vectors to operate on inside dt

# define these values we're interested in using lapply and .SD/.SDcols
# .SD refers to the columns that we're interested in applying the function to
# note that we could have also applied the paste function directly inside the dt call
# here, we are adding a fourth argument to the inside of []
rdt[, (max_cols) := lapply(.SD, function(x) max(x, na.rm = TRUE)), .SDcols=cols, by = subreddit]

# can use setorder to see the result by subreddit
setorder(rdt, subreddit)

# can set orders by multiple columns as well
setorderv(rdt, c("subreddit", "username"))

# can also set column order, but need to supply all columns
setcolorder(rdt, rev(colnames(rdt)))

# we can use rbindlist to rbind two data tables, ignoring the column names of the second and just using the first
# useful when we want to compare data objects with different names
rbindlist(list(red, setcolorder(rdt, rev(colnames(rdt)))))
# rbindlist takes a list of data tables or lists as arguments

# special operators that are useful with by
# .SD refers to a data.table that is the subset of the data for the by group (excluding the columns in by)
# .BY refers to a list with a vector of length 1 with each item in by (can be helpful to return if generating "by" dynamiicaly)
# .N is the number of rows in each by group
# .I is the row position of each element in the by group
# .GRP is a numerical identifier for each group

# assignemnt with list -- a way to remove columns
# when we add columns, we use :=
# if we only want specific columns back, we can use list
rdt_uservars = rdt[, list(username = username, subreddit = subreddit, count = count, min = min, avg_comm = avg_comm, var_comm = var_comm, skew_comm = skew_comm)]

# we can combine list with unique to collapse and aggregate the data
# since these summary stat columns are the same within users, we can write
rdt_uservars = unique(rdt[, list(username = username, avg_comm = avg_comm, var_comm = var_comm, skew_comm = skew_comm)])
# to get these measures in a table at the user level
nrow(rdt_uservars) == length(unique(rdt_uservars$username))

# take a few minutes to try to make the following and discuss
# let's get a subreddit-level measure of mean commenting activity on each subreddit by each user
# but let's only generate these measures for users that actually comment on the subreddit
# to get a sense of how active the members that really participate are
# a conditional mean similar to the mean/median tie strength calculations from exercise 1

# exercise to try on your own
# the utc time variable is given in epoch time, or the number of seconds that have elapsed since 12am, January 1st, 1970
# create a data table that illustrates the mean, variance, and skew of each users' commenting behavior for each month that they are active
# are a current month's summary statistics strongly related to the summary statistics from the previous month?

