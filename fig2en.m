function fig2en
h = gcf;
saveas(h,'enfigtmp.fig');
saveas(h,'enfigtmp.png');
% Define these variables appropriately:
myaddress = 'lab2qo@gmail.com';
mypassword = 'cesium133';
% preferences
setpref('Internet','E_mail',myaddress);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',myaddress);
setpref('Internet','SMTP_Password',mypassword);
% crazy java properties
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', ...
    'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');
attfile = {'enfigtmp.fig','enfigtmp.png'};
sendmail('acmcclung.8c279@m.evernote.com','new figure @M3P','',attfile);
delete('enfigtmp.fig','enfigtmp.png');
end
