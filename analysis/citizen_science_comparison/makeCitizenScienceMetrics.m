function makeCitizenScienceMetrics(datapath,outputdir,outfile,ccp_tutorialdata,mmos_stats)
%These metrics are taken from:
%Defining and Measuring Success in Online Citizen Science: A Case Study of Zooniverse Projects
%http://eprints.whiterose.ac.uk/86535/1/Success%20in%20Citizen%20Science%20Final%20Submission_Revised%20Version%20Feb%202015.pdf
%
%Written by: Devin Sullivan 

if nargin<1
    datapath = './parsed_data_201607121468334823219_tasks.mat';
end

if nargin<2 || isempty(outputdir)
    outputdir = '../results';
end

%%Defined parameters 
if nargin<3 || isempty(outfile)
    outfile = 'PD_citizenscience_metrics.txt';
end

if nargin<4 || isempty(ccp_tutorialdata)
    ccp_tutorialdata = '../data/ccp_data/tutorialstartfinish.xlsx';
end
if nargin<5 || isempty(mmos_stats)
    mmos_stats = '../data/mmos_data/mmos_stats.txt';
end
%constants
ndays_year = 365;
full_time_hrs = 35;%hours per week defined by cited paper
num_citations = 0;

[number_of_seconds,tot_players,...
    mean_per_player,median_per_player,...
    players_with_no_tasks] = getMMOSStats(mmos_stats);
% number_of_seconds = 420543320;%**not counting accounts open for more than 2 minutes

number_of_hours = number_of_seconds/60/60;
%players_tutorial = NaN;%number of players that played the tutorial only
%number of players total (tutorial + those past it)
%tot_players = 164783;%This is the number of capsuleers, so true players is less
% tot_players_toAug = 328181;
% tot_through_tut_toAug = 72263;
% players_tutorial_toAug = tot_players_toAug-tot_through_tut_toAug;
papers_w_cscientists = 0;
num_published = 2; %see list below 

% mean_per_player = 106.98;
% median_per_player = 9;
mean_active_period = (number_of_hours/tot_players)/ndays_year;%10/ndays_year;
mean_vote_per_task = 13;

%Cell Atlas Paper - in review
%Project Discovery Paper - in preparation

%twitter & instagram data taken from keyhole.co
%@HPA_Discovery - 265 (October 28, 2016)
hpa_discovery = 265;
humanproteome = 962;
humanproteinatlas_insta = 33;%instagram
humanproteinatlas_fb = 20;%rough estimate
devinpsullivan = 343;
%manually counted HPA IOTW blogs 
hpa_blog = 11;
%Reddit stats taken from http://www.redective.com/ (October 28, 2016)
hpa_dichroic_eve_submissions = 4;
hpa_dichroic_eve_comments = 163;
hpa_dichroic_pd_submissions = 27;
hpa_dichroic_pd_comments = 109;
hpa_illuminator_eve_submissions = 1;
hpa_illuminator_eve_comments = 328;
hpa_illuminator_pd_submissions = 7;
hpa_illuminator_pd_comments = 167;
pd_submissions_tot = 162;%manually acquired from /r/projectdiscovery
total_redditHPA_submissions = hpa_dichroic_pd_submissions+...
    hpa_dichroic_eve_submissions+hpa_illuminator_eve_submissions+...
    hpa_illuminator_pd_submissions;
total_reddit_submissions = pd_submissions_tot+...
    hpa_dichroic_eve_submissions+hpa_illuminator_eve_submissions;
total_redditHPA_comments = hpa_dichroic_eve_comments+...
    hpa_dichroic_pd_comments+hpa_illuminator_eve_comments+hpa_illuminator_pd_comments;
%%youtube 
hpayoutube_account = 9;
youtube_view = 62+57+204+1560+131+39786+2771+8626+4511;

%EVE forum posts 
%%%CAN'T DO THIS ONE RIGHT NOW, ESTIMATING%%%
forum_posts = 3;
forum_comments_hpa = 15;
forum_comments_tot = 250;
%%%
%Player generated content 
%%%CAN'T DO THIS RIGHT NOW, ESTIMATING%%%
podcasts = 4;
%%%
%%Official talks posted 
talks = 10;%rough estimate

%General metrics 
load(datapath,'runningdates')
num_days = max(runningdates);
project_age = (num_days/ndays_year);

%Contribution to Science 

  %Data Value 
  
  publication_rate = num_published/project_age^2;
%   completeness = length(originalCode)/length(unique(originalCode));
  completeness = length(originalCode)/length(unique(originalCode))*mean_vote_per_task;
  academic_impact = num_citations/project_age^2;%undefined still
  
  
  %Project design 
  
%   resource_saving = 1-((num_days*24)/(number_of_hours/35));
  resource_saving = 1-((num_days/7)/(number_of_hours/full_time_hrs));
  distribution_effort = NaN;%don't know what the Gini coeff is...
  effective_training = computeEffectiveTraining(ccp_tutorialdata);
  %1-0.47;%this 0.47 is the mean computed from Hjalti %(players_tutorial/tot_players);
  
  
%Public engagement 
  
  %Dissemination 
  collaboration = papers_w_cscientists/project_age^2;
  
  tot_communication = sum([hpa_discovery,humanproteome,...
      humanproteinatlas_insta,humanproteinatlas_fb,devinpsullivan,...
      hpa_blog,hpa_dichroic_eve_submissions,hpa_dichroic_eve_comments,...
      hpa_dichroic_pd_submissions,hpa_dichroic_pd_comments,...
      hpa_illuminator_eve_submissions,hpa_illuminator_eve_comments,...
      hpa_illuminator_pd_submissions,hpa_illuminator_pd_comments,...
      hpayoutube_account,forum_posts,forum_comments_hpa,...
      podcasts,talks]);
  communication = tot_communication/project_age^2;
  tot_interaction = sum([hpa_dichroic_eve_comments,...
      hpa_dichroic_pd_comments,hpa_illuminator_eve_comments,...
      hpa_illuminator_pd_comments,forum_comments_hpa,...
      podcasts,talks]);
  interaction = tot_interaction/project_age^2;
  
  
  %Participation 
  project_appeal = tot_players/project_age^2;
  %using mean instead of median because I assume we can't get median
  sustained_engagement = mean_active_period/project_age^2; 
  public_contribution = mean_per_player/project_age^2;
  
  
  %%%end
  metricnames = {'publication rate','completeness of analysis','data impact',...
      'resource savings','distribution of effort','effective training',...
      'collaboration','communication','interaction','project appeal',...
      'sustained engagement','public contribution'};
  metrics = [publication_rate,completeness,academic_impact,resource_saving,...
      distribution_effort,effective_training,collaboration,communication,...
      interaction,project_appeal,sustained_engagement,public_contribution];
  
  %write out results 
  fid = fopen([outputdir,filesep,outfile],'w')
  headers = {'measurement','proxy'};
  fprintf(fid,'%s,%s\n',headers{:});
  for i = 1:length(metricnames)
      
      fprintf(fid,'"%s",%f\n',metricnames{i},metrics(i));
  end
  fclose(fid)   